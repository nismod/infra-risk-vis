-- irv-metrics: relational data model (PostgreSQL 14+ with PostGIS 3).
--
-- Applied by the loaders (irv_metrics.loaders.common.apply_schema) when the
-- `metrics` schema is missing. Rationale: docs/data-model/README.md,
-- docs/integration-plan.md ("Boundaries", "Data cube storage"); public contract:
-- docs/api-v1.md; mapping from the upstream GRI API: backend/irv_metrics/specs.py.
--
-- Shape: a statistical data cube.
--   reference   country, unit, section
--   boundaries  boundary_set, boundary_set_level, area (geometry), area_country
--   vocabulary  dimension, dimension_member, member_rollup
--   structure   dataset, dataset_alias, metric, dataset_dimension, dataset_dimension_member
--   facts       series (dataset + metric + one member per dimension), observation, curve_point
--   provenance  load_record (one record per dataset x boundary set x country)
--
-- Every API response kind (series, profile, ladder) is a query over observation:
--   series  = one series, every area of a set at one level related to a country
--   profile = one area, several series that differ in one "breakdown" dimension
--   ladder  = one area, series that differ only in return_period_years
--
-- Not partitioned (see docs/integration-plan.md, "Recommendation").

CREATE EXTENSION IF NOT EXISTS postgis;
CREATE SCHEMA IF NOT EXISTS metrics;
SET search_path = metrics, public;

-- ---------------------------------------------------------------------------
-- Reference data
-- ---------------------------------------------------------------------------

CREATE TABLE unit (
    unit_id     text PRIMARY KEY,                        -- API unit string: km2, usd_per_year, ...
    label       text NOT NULL,                           -- display: km², USD/year
    quantity    text NOT NULL                            -- area, length, count, currency, ratio, ...
);

CREATE TABLE section (
    section_id  text PRIMARY KEY,                        -- hazard, exposure, ... (API family)
    title       text NOT NULL,
    description text,
    sort_order  smallint NOT NULL
);

-- ---------------------------------------------------------------------------
-- Boundaries (docs/api-v1.md, "Boundary sets and areas")
-- ---------------------------------------------------------------------------

-- One row per geometry release or grouping, e.g. cgaz-6.0.0, hydrobasins-v1c,
-- wb-income-fy26. kind says how its areas relate to countries:
--   admin    nested, partitions each country; level 0 is the country (area_id = ISO3)
--   overlay  drawn without regard to borders; related to every country it overlaps
--            by at least min_share of its area, by intersection with level 0 of
--            the admin set countries_from (computed at load)
--   group    named sets of countries, no geometry; related to each member
CREATE TABLE boundary_set (
    boundary_set    text PRIMARY KEY CHECK (boundary_set ~ '^[a-z][a-z0-9.-]*$'),
    kind            text NOT NULL CHECK (kind IN ('admin', 'overlay', 'group')),
    title           text NOT NULL,
    attribution     text,
    attribution_url text,
    countries_from  text REFERENCES boundary_set,         -- overlay only: an admin set
    min_share       double precision CHECK (min_share BETWEEN 0 AND 1),  -- overlay only
    CHECK ((kind = 'overlay') = (countries_from IS NOT NULL)),
    CHECK ((kind = 'overlay') OR min_share IS NULL)
);

-- Levels belong to the set: admin sets number from 0 (the country); overlay and
-- group sets keep their own numbers (HydroBASINS 1-12).
CREATE TABLE boundary_set_level (
    boundary_set text NOT NULL REFERENCES boundary_set ON DELETE CASCADE,
    level        smallint NOT NULL CHECK (level >= 0),
    label        text,                                   -- local name: County, Parish, Sub-basin
    PRIMARY KEY (boundary_set, level)
);

CREATE TABLE country (
    iso3         char(3) PRIMARY KEY CHECK (iso3 = upper(iso3)),
    name         text NOT NULL,
    boundary_set text REFERENCES boundary_set            -- the country's admin set (map and boundaries default)
);

-- Every area of every set. area_id is the published identifier, unique within
-- its set across levels: ISO3 at admin level 0, the CGAZ shapeID at levels 1-2
-- (e.g. 32016919B10079274489584), a Pfafstetter code for a basin. area_key is
-- the internal integer that observations are keyed on.
CREATE TABLE area (
    area_key        int4 GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
    boundary_set    text NOT NULL REFERENCES boundary_set,
    area_id         text NOT NULL,
    level           smallint NOT NULL CHECK (level >= 0),
    name            text NOT NULL,
    parent_area_key int4 REFERENCES area,                -- same set; derived by spatial join on load when not given
    geometry        geometry(MultiPolygon, 4326),        -- NULL for group sets
    area_km2        double precision GENERATED ALWAYS AS (ST_Area(geometry::geography) / 1e6) STORED,
    source          text,                                -- e.g. 'GRI API /boundaries/ken/1'
    UNIQUE (boundary_set, area_id)
);
CREATE INDEX area_set_level ON area (boundary_set, level);
CREATE INDEX area_parent ON area (parent_area_key);
CREATE INDEX area_geometry ON area USING gist (geometry);

-- Which countries each area is related to: the country scope of every
-- /countries/{iso3}/... route is a filter through this table.
--   admin    one row, share 1
--   overlay  one row per overlapped country, share = fraction of the area in it
--   group    one row per member country, share 1
CREATE TABLE area_country (
    area_key int4 NOT NULL REFERENCES area ON DELETE CASCADE,
    iso3     char(3) NOT NULL REFERENCES country,
    share    double precision NOT NULL CHECK (share > 0 AND share <= 1),
    PRIMARY KEY (area_key, iso3)
);
CREATE INDEX area_country_iso3 ON area_country (iso3);

-- ---------------------------------------------------------------------------
-- Controlled vocabulary: one row per dimension and per member. The member_id is
-- the canonical slug (the API value where the API has one). Labels live here,
-- once; per-dataset wording differences go in dataset_dimension_member.
-- ---------------------------------------------------------------------------

CREATE TABLE dimension (
    dimension_id text PRIMARY KEY,                       -- demographic_group, return_period_years, ...
    label        text NOT NULL,
    description  text,
    is_ordered   boolean NOT NULL DEFAULT false
);

CREATE TABLE dimension_member (
    dimension_id  text NOT NULL REFERENCES dimension,
    member_id     text NOT NULL,                         -- q1, childbearing, 100, cat2plus, ...
    label         text NOT NULL,                         -- 'Female Aged 15–49'
    short_label   text,                                  -- 'Females 15–49' (chart axes, legends)
    description   text,
    sort_order    smallint NOT NULL,
    numeric_value double precision,                      -- 100 for RP100, 15 for pop_within_15min
    code          text,                                  -- external codes, e.g. GHS-SMOD '30'
    PRIMARY KEY (dimension_id, member_id)
);

-- Aggregate members and how they are built: total = sum of sectors,
-- bottom_40 = q1 + q2, cat2plus = cat2 + ... + cat5. Lets analytics recompute or
-- validate aggregates rather than trusting them blindly.
CREATE TABLE member_rollup (
    dimension_id text NOT NULL,
    parent_id    text NOT NULL,
    child_id     text NOT NULL,
    PRIMARY KEY (dimension_id, parent_id, child_id),
    FOREIGN KEY (dimension_id, parent_id) REFERENCES dimension_member,
    FOREIGN KEY (dimension_id, child_id) REFERENCES dimension_member
);

-- ---------------------------------------------------------------------------
-- Dataset structure ("data structure definition")
-- ---------------------------------------------------------------------------

-- One dataset per map layer / indicator card. dataset_id is the API path
-- fragment 'family/layer' (hazard/river-flooding, risk/social-accessibility).
CREATE TABLE dataset (
    dataset_id         text PRIMARY KEY,
    section_id         text NOT NULL REFERENCES section,
    title              text NOT NULL,
    description        text,
    documentation_path text,                             -- /docs/hazard/river-flooding
    source_note        text,                             -- e.g. 'JRC return-period flood maps'
    status             text NOT NULL DEFAULT 'available' CHECK (status IN ('available', 'unavailable')),
    sort_order         smallint NOT NULL
);

-- The three identifier systems in play (see README "Identifiers").
CREATE TABLE dataset_alias (
    alias       text NOT NULL,
    system      text NOT NULL CHECK (system IN ('bundle_layer', 'reconstruction_card')),
    dataset_id  text NOT NULL REFERENCES dataset,
    PRIMARY KEY (system, alias)
);

-- Measures. aggregation says how values combine across areas (or not at all):
-- counts, lengths, people, USD sum; depths/percentages/ratios do not.
CREATE TABLE metric (
    dataset_id   text NOT NULL REFERENCES dataset,
    metric_id    text NOT NULL,                          -- flooded_area, average_annual_loss, ...
    label        text NOT NULL,
    unit_id      text NOT NULL REFERENCES unit,
    aggregation  text NOT NULL CHECK (aggregation IN ('sum', 'mean', 'weighted_mean', 'none')),
    weight_metric text,                                  -- for weighted_mean, e.g. population
    description  text,
    PRIMARY KEY (dataset_id, metric_id)
);

-- Which dimensions a dataset's series are keyed by.
CREATE TABLE dataset_dimension (
    dataset_id     text NOT NULL REFERENCES dataset,
    dimension_id   text NOT NULL REFERENCES dimension,
    default_member text,
    -- When the dimension applies, e.g. {"metric": ["per_100k"]} for demographic_group
    -- on facilities, or {"infrastructure_type": ["road"]} for road_class on direct damage.
    applies_when   jsonb,
    PRIMARY KEY (dataset_id, dimension_id),
    FOREIGN KEY (dimension_id, default_member) REFERENCES dimension_member
);

-- The members valid for a dataset (e.g. wealth_group includes bottom_40 for
-- vulnerability but total for adaptation outcomes), with optional wording overrides
-- ('walking' is 'Pedestrian'; 'total' is 'All' for road class).
CREATE TABLE dataset_dimension_member (
    dataset_id     text NOT NULL,
    dimension_id   text NOT NULL,
    member_id      text NOT NULL,
    label_override text,
    PRIMARY KEY (dataset_id, dimension_id, member_id),
    FOREIGN KEY (dataset_id, dimension_id) REFERENCES dataset_dimension,
    FOREIGN KEY (dimension_id, member_id) REFERENCES dimension_member
);

-- ---------------------------------------------------------------------------
-- Facts
-- ---------------------------------------------------------------------------

-- A series is one fully specified slice: dataset + metric + one member for each
-- applicable dimension. dims is the canonical {dimension: member} map. series_member
-- holds the same pairs with foreign keys, so every member is a validated vocabulary
-- entry allowed for the dataset. The loader writes both; a check query compares them.
-- A series is not tied to a boundary set: the same series may have values on
-- several sets (counties, countries and basins).
CREATE TABLE series (
    series_id  int4 GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
    dataset_id text NOT NULL,
    metric_id  text NOT NULL,
    dims       jsonb NOT NULL DEFAULT '{}'::jsonb CHECK (jsonb_typeof(dims) = 'object'),
    FOREIGN KEY (dataset_id, metric_id) REFERENCES metric,
    UNIQUE (dataset_id, metric_id, dims),
    UNIQUE (series_id, dataset_id)                       -- target for series_member's composite FK
);
CREATE INDEX series_dims ON series USING gin (dims jsonb_path_ops);

CREATE TABLE series_member (
    series_id    int4 NOT NULL,
    dataset_id   text NOT NULL,
    dimension_id text NOT NULL,
    member_id    text NOT NULL,
    PRIMARY KEY (series_id, dimension_id),
    FOREIGN KEY (series_id, dataset_id) REFERENCES series (series_id, dataset_id) ON DELETE CASCADE,
    FOREIGN KEY (dataset_id, dimension_id, member_id) REFERENCES dataset_dimension_member
);

-- One value per series per area: the slim long table. No country column (the
-- country is reached through area_country) and no per-row load id (see load_record).
-- Values are stored as published for each level: national values are NOT sums
-- of regional values in the source data (zonal statistics are computed
-- independently per level; see docs/data-model/README.md).
CREATE TABLE observation (
    series_id int4 NOT NULL REFERENCES series ON DELETE CASCADE,
    area_key  int4 NOT NULL REFERENCES area,
    value     double precision,                          -- NULL = published as null
    PRIMARY KEY (series_id, area_key)
);
CREATE INDEX observation_area ON observation (area_key);

-- Concentration curves (Lorenz-type: cumulative exposure share vs cumulative
-- population share ranked by wealth). National only (area = admin level 0).
CREATE TABLE curve_point (
    series_id int4 NOT NULL REFERENCES series ON DELETE CASCADE,
    area_key  int4 NOT NULL REFERENCES area,
    seq       smallint NOT NULL,
    x         double precision NOT NULL,                 -- cumulative_population_share
    y         double precision NOT NULL,                 -- cumulative exposure share
    PRIMARY KEY (series_id, area_key, seq)
);

-- Provenance: one record per (dataset, boundary set, country), replaced on
-- reload. iso3 is NULL for a set loaded globally (an overlay or group set, or a
-- global admin file). The rows a record covers are the dataset's observations
-- on the set's areas related to that country (all of the set's areas when iso3
-- is NULL). finished_at changes whenever a load completes; the API derives its
-- cache version and ETags from it.
CREATE TABLE load_record (
    load_id      int4 GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
    dataset_id   text NOT NULL REFERENCES dataset,
    boundary_set text NOT NULL REFERENCES boundary_set,
    iso3         char(3) REFERENCES country,
    source       text NOT NULL,                          -- 'gri-api scrape KEN', 'csv api-data/KEN_adm1_...'
    started_at   timestamptz NOT NULL DEFAULT now(),
    finished_at  timestamptz,
    notes        jsonb
);
-- One record per key, NULL iso3 included (UNIQUE NULLS NOT DISTINCT needs PG 15).
CREATE UNIQUE INDEX load_record_key ON load_record (dataset_id, boundary_set, coalesce(iso3, ''));
CREATE INDEX load_record_country ON load_record (iso3);

-- ---------------------------------------------------------------------------
-- Views and query shapes behind the public API
-- ---------------------------------------------------------------------------

-- Long, analysis-friendly view: one row per value and related country. Overlay
-- and group areas appear once per country they relate to (with its share).
CREATE VIEW observation_long AS
SELECT s.dataset_id, s.metric_id, s.dims, a.boundary_set, ac.iso3, ac.share, a.level, a.area_id, a.name,
       o.value, m.unit_id
FROM observation o
JOIN series s USING (series_id)
JOIN metric m USING (dataset_id, metric_id)
JOIN area a USING (area_key)
JOIN area_country ac USING (area_key);

-- series:  GET /metrics/hazard/river-flooding/ken/1?metric=flooded_area&return_period_years=10
--   SELECT area_id, name, value, unit_id FROM observation_long
--   WHERE dataset_id = 'hazard/river-flooding' AND metric_id = 'flooded_area'
--     AND dims = '{"return_period_years": "10"}' AND boundary_set = 'cgaz-6.0.0'
--     AND iso3 = 'KEN' AND level = 1;
--
-- profile: GET /metrics/hazard/river-flooding/ken/profile/{area_id}?return_period_years=10
--   SELECT metric_id, value, unit_id FROM observation_long
--   WHERE dataset_id = 'hazard/river-flooding' AND boundary_set = 'cgaz-6.0.0'
--     AND area_id = :area_id AND iso3 = 'KEN' AND dims = '{"return_period_years": "10"}';
--
-- ladder:  GET /metrics/hazard/river-flooding/ken/ladder/{area_id}?metric=flooded_area
--   SELECT dims->>'return_period_years' AS rp, value FROM observation_long
--   WHERE dataset_id = 'hazard/river-flooding' AND metric_id = 'flooded_area'
--     AND boundary_set = 'cgaz-6.0.0' AND area_id = :area_id AND iso3 = 'KEN'
--     AND dims - 'return_period_years' = '{}'
--   ORDER BY (dims->>'return_period_years')::int;
--
-- display_mode=percentage for population/wealth is derived, not stored:
--   100 * value / value_of(demographic_group = 'total'), rounded to 3 decimals
--   (verified against the scrape).
