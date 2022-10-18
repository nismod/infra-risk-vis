import { colorMap } from '@/lib/color-map';
import { StyleParams, ViewLayer } from '@/lib/data-map/view-layers';
import { dataColorMap } from '@/lib/deck/props/color-map';
import { border, lineStyle, pointRadius } from '@/lib/deck/props/style';
import { fillColor, strokeColor } from '@/lib/deck/props/style';
import { makeConfig } from '@/lib/helpers';

import { getAssetDataAccessor } from '@/config/assets/data-access';
import { COLORS } from '@/config/colors';

import { infrastructureViewLayer } from './infrastructure-view-layer';

function infraStyle(layer: string, defaultStyle, styleParams: StyleParams) {
  if (styleParams?.colorMap) {
    const { fieldSpec, colorSpec } = styleParams.colorMap;
    return dataColorMap(getAssetDataAccessor(layer, fieldSpec), colorMap(colorSpec));
  }
  return defaultStyle;
}

enum RoadClass {
  class_a = 'class_a',
  class_b = 'class_b',
  class_c = 'class_c',
  metro = 'metro',
  other = 'other',
  track = 'track',
}

const roadClassLookup = {
  road_edges_class_a: RoadClass.class_a,
  road_edges_class_b: RoadClass.class_b,
  road_edges_class_c: RoadClass.class_c,
  road_edges_metro: RoadClass.metro,
  road_edges_track: RoadClass.track,
  road_edges_other: RoadClass.other,
};

const roadColor = {
  [RoadClass.class_a]: COLORS.roads_class_a.css,
  [RoadClass.class_b]: COLORS.roads_class_b.css,
  [RoadClass.class_c]: COLORS.roads_class_c.css,
  [RoadClass.metro]: COLORS.roads_metro.css,
  [RoadClass.track]: COLORS.roads_unknown.css,
  [RoadClass.other]: COLORS.roads_unknown.css,
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
  infrastructureViewLayer('elec_edges_high', ({ zoom, styleParams }) => [
    strokeColor(infraStyle('elec_edges_high', COLORS.electricity_high.deck, styleParams)),
    lineStyle(zoom),
  ]),
  infrastructureViewLayer('elec_edges_low', ({ zoom, styleParams }) => [
    strokeColor(infraStyle('elec_edges_low', COLORS.electricity_low.deck, styleParams)),
    lineStyle(zoom),
  ]),
  electricitySourceViewLayer('elec_nodes_diesel'),
  electricitySourceViewLayer('elec_nodes_gas'),
  electricitySourceViewLayer('elec_nodes_hydro'),
  electricitySourceViewLayer('elec_nodes_solar'),
  electricitySourceViewLayer('elec_nodes_wind'),
  infrastructureViewLayer('elec_nodes_demand', ({ zoom, styleParams }) => [
    border(),
    fillColor(infraStyle('elec_nodes_demand', COLORS.electricity_low.deck, styleParams)),
    pointRadius(zoom),
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
  roadsViewLayer('road_edges_class_a'),
  roadsViewLayer('road_edges_class_b'),
  roadsViewLayer('road_edges_class_c'),
  roadsViewLayer('road_edges_metro'),
  roadsViewLayer('road_edges_track'),
  roadsViewLayer('road_edges_other'),
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
