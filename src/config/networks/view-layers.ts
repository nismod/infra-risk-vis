import { makeConfig } from 'lib/helpers';
import { border, lineStyle, pointRadius } from 'lib/deck/props/style';

import { COLORS } from '../colors';
import { infrastructureViewLayer } from './infrastructure-view-layer';
import { StyleParams, ViewLayer } from 'lib/data-map/view-layers';
import { fillColor, strokeColor } from 'lib/deck/props/style';
import { dataColorMap } from 'lib/deck/props/color-map';
import { getAssetDataAccessor } from 'config/assets/data-access';
import { colorMap } from 'lib/color-map';

function infraStyle(layer: string, defaultStyle, styleParams: StyleParams) {
  if (styleParams?.colorMap) {
    const { fieldSpec, colorSpec } = styleParams.colorMap;
    return dataColorMap(getAssetDataAccessor(layer, fieldSpec), colorMap(colorSpec));
  }
  return defaultStyle;
}

enum RoadClass {
  motorway = 'motorway',
  trunk = 'trunk',
  primary = 'primary',
  secondary = 'secondary',
  tertiary = 'tertiary',
}

const roadClassLookup = {
  road_edges_motorway: RoadClass.motorway,
  road_edges_trunk: RoadClass.trunk,
  road_edges_primary: RoadClass.primary,
  road_edges_secondary: RoadClass.secondary,
  road_edges_tertiary: RoadClass.tertiary,
};

const roadColor = {
  [RoadClass.motorway]: COLORS.roads_motorway.css,
  [RoadClass.trunk]: COLORS.roads_trunk.css,
  [RoadClass.primary]: COLORS.roads_primary.css,
  [RoadClass.secondary]: COLORS.roads_secondary.css,
  [RoadClass.tertiary]: COLORS.roads_tertiary.css,
};
function roadsViewLayer(asset_id) {
  return infrastructureViewLayer(asset_id, ({ zoom, styleParams }) => [
    strokeColor(
      infraStyle(
        asset_id,
        dataColorMap(
          () => asset_id,
          (x) => roadColor[roadClassLookup[x]],
        ),
        styleParams,
      ),
    ),
    lineStyle(zoom),
  ]);
}

function potableNodesViewLayer(asset_id) {
  return infrastructureViewLayer(asset_id, ({ zoom, styleParams }) => [
    border(),
    pointRadius(zoom),
    fillColor(infraStyle(asset_id, COLORS.water_supply.deck, styleParams)),
  ]);
}

function wastewaterNodesViewLayer(asset_id) {
  return infrastructureViewLayer(asset_id, ({ zoom, styleParams }) => [
    border(),
    pointRadius(zoom),
    fillColor(infraStyle(asset_id, COLORS.water_wastewater.deck, styleParams)),
  ]);
}

function electricitySourceViewLayer(asset_id) {
  return infrastructureViewLayer(asset_id, ({ zoom, styleParams }) => [
    border(),
    fillColor(infraStyle(asset_id, COLORS.electricity_high.deck, styleParams)),
    pointRadius(zoom),
  ]);
}

export const INFRASTRUCTURE_VIEW_LAYERS = makeConfig<ViewLayer, string>([
  infrastructureViewLayer('elec_edges', ({ zoom, styleParams }) => [
    strokeColor(infraStyle('elec_edges', COLORS.electricity_high.deck, styleParams)),
    lineStyle(zoom),
  ]),
  infrastructureViewLayer('elec_edges_high', ({ zoom, styleParams }) => [
    strokeColor(infraStyle('elec_edges_high', COLORS.electricity_high.deck, styleParams)),
    lineStyle(zoom),
  ]),
  infrastructureViewLayer('elec_edges_low', ({ zoom, styleParams }) => [
    strokeColor(infraStyle('elec_edges_low', COLORS.electricity_low.deck, styleParams)),
    lineStyle(zoom),
  ]),
  electricitySourceViewLayer('elec_nodes_biomass'),
  electricitySourceViewLayer('elec_nodes_coal'),
  electricitySourceViewLayer('elec_nodes_diesel'),
  electricitySourceViewLayer('elec_nodes_gas'),
  electricitySourceViewLayer('elec_nodes_getertiarymal'),
  electricitySourceViewLayer('elec_nodes_hydro'),
  electricitySourceViewLayer('elec_nodes_oil'),
  electricitySourceViewLayer('elec_nodes_solar'),
  electricitySourceViewLayer('elec_nodes_wind'),
  infrastructureViewLayer('elec_areas_demand', ({ zoom, styleParams }) => [
    border(),
    fillColor(infraStyle('elec_nodes_demand', COLORS.electricity_unknown.deck, styleParams))
  ]),
  infrastructureViewLayer('elec_nodes_pole', ({ zoom, styleParams }) => [
    border(),
    fillColor(infraStyle('elec_nodes_pole', COLORS.electricity_unknown.deck, styleParams)),
    pointRadius(zoom),
  ]),
  infrastructureViewLayer('elec_nodes_substation', ({ zoom, styleParams }) => [
    border(),
    fillColor(infraStyle('elec_nodes_substation', COLORS.electricity_unknown.deck, styleParams)),
    pointRadius(zoom),
  ]),
  infrastructureViewLayer('rail_edges', ({ zoom, styleParams }) => [
    strokeColor(infraStyle('rail_edges', COLORS.railway.deck, styleParams)),
    lineStyle(zoom),
  ]),
  infrastructureViewLayer('rail_stations', ({ zoom, styleParams }) => [
    border(),
    fillColor(infraStyle('rail_stations', COLORS.railway.deck, styleParams)),
    pointRadius(zoom),
  ]),
  infrastructureViewLayer('rail_junctions', ({ zoom, styleParams }) => [
    border(),
    fillColor(infraStyle('rail_junctions', COLORS.railway.deck, styleParams)),
    pointRadius(zoom),
  ]),
  roadsViewLayer('road_edges_motorway'),
  roadsViewLayer('road_edges_trunk'),
  roadsViewLayer('road_edges_primary'),
  roadsViewLayer('road_edges_secondary'),
  roadsViewLayer('road_edges_tertiary'),
  infrastructureViewLayer('road_bridges', ({ zoom, styleParams }) => [
    border(),
    fillColor(infraStyle('road_bridges', COLORS.bridges.deck, styleParams)),
    pointRadius(zoom),
  ]),
  infrastructureViewLayer('airport_runways', ({ zoom, styleParams }) => [
    border([0, 0, 0]),
    fillColor(infraStyle('airport_runways', COLORS.airports.deck, styleParams)),
  ]),
  infrastructureViewLayer('airport_terminals', ({ zoom, styleParams }) => [
    border([0, 0, 0]),
    fillColor(infraStyle('airport_terminals', COLORS.airports.deck, styleParams)),
  ]),
  infrastructureViewLayer('port_areas_break', ({ zoom, styleParams }) => [
    border(),
    fillColor(infraStyle('port_areas_break', COLORS.ports.deck, styleParams)),
  ]),
  infrastructureViewLayer('port_areas_container', ({ zoom, styleParams }) => [
    border(),
    fillColor(infraStyle('port_areas_container', COLORS.ports.deck, styleParams)),
  ]),
  infrastructureViewLayer('port_areas_industry', ({ zoom, styleParams }) => [
    border(),
    fillColor(infraStyle('port_areas_industry', COLORS.ports.deck, styleParams)),
  ]),
  infrastructureViewLayer('port_areas_silo', ({ zoom, styleParams }) => [
    border(),
    fillColor(infraStyle('port_areas_silo', COLORS.ports.deck, styleParams)),
  ]),
  potableNodesViewLayer('water_potable_nodes_booster'),
  potableNodesViewLayer('water_potable_nodes_catchment'),
  potableNodesViewLayer('water_potable_nodes_entombment'),
  potableNodesViewLayer('water_potable_nodes_filter'),
  potableNodesViewLayer('water_potable_nodes_intake'),
  potableNodesViewLayer('water_potable_nodes_well'),
  potableNodesViewLayer('water_potable_nodes_pump'),
  potableNodesViewLayer('water_potable_nodes_relift'),
  potableNodesViewLayer('water_potable_nodes_reservoir'),
  potableNodesViewLayer('water_potable_nodes_river_source'),
  potableNodesViewLayer('water_potable_nodes_spring'),
  potableNodesViewLayer('water_potable_nodes_tank'),
  potableNodesViewLayer('water_potable_nodes_sump'),
  potableNodesViewLayer('water_potable_nodes_tp'),

  infrastructureViewLayer('water_potable_edges', ({ zoom, styleParams }) => [
    lineStyle(zoom),
    strokeColor(infraStyle('water_potable_edges', COLORS.water_supply.deck, styleParams)),
  ]),
  infrastructureViewLayer('water_irrigation_edges', ({ zoom, styleParams }) => [
    lineStyle(zoom),
    strokeColor(infraStyle('water_irrigation_edges', COLORS.water_irrigation.deck, styleParams)),
  ]),
  infrastructureViewLayer('water_irrigation_nodes', ({ zoom, styleParams }) => [
    border(),
    pointRadius(zoom),
    fillColor(infraStyle('water_irrigation_nodes', COLORS.water_irrigation.deck, styleParams)),
  ]),
  infrastructureViewLayer('water_waste_sewer_gravity', ({ zoom, styleParams }) => [
    lineStyle(zoom),
    strokeColor(infraStyle('water_waste_sewer_gravity', COLORS.water_wastewater.deck, styleParams)),
  ]),
  infrastructureViewLayer('water_waste_sewer_pressure', ({ zoom, styleParams }) => [
    lineStyle(zoom),
    strokeColor(infraStyle('water_waste_sewer_pressure', COLORS.water_wastewater.deck, styleParams)),
  ]),
  wastewaterNodesViewLayer('water_waste_nodes_sump'),
  wastewaterNodesViewLayer('water_waste_nodes_pump'),
  wastewaterNodesViewLayer('water_waste_nodes_relift'),
  wastewaterNodesViewLayer('water_waste_nodes_wwtp'),
]);
