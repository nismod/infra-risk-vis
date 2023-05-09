import { makeConfig } from 'lib/helpers';
import { ScaleLevel, border, iconColor, iconSize, lineStyle, pointRadius } from 'lib/deck/props/style';

import { COLORS } from '../colors';
import { infrastructureViewLayer } from './infrastructure-view-layer';
import { StyleParams, ViewLayer } from 'lib/data-map/view-layers';
import { fillColor, strokeColor } from 'lib/deck/props/style';
import { dataColorMap } from 'lib/deck/props/color-map';
import { getAssetDataAccessor } from 'config/assets/data-access';
import { colorMap } from 'lib/color-map';
import { NETWORKS_METADATA } from './metadata';
import { iconType } from 'lib/map-shapes/deck-icon';

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
  [RoadClass.class_a]: COLORS.roads_class_a.deck,
  [RoadClass.class_b]: COLORS.roads_class_b.deck,
  [RoadClass.class_c]: COLORS.roads_class_c.deck,
  [RoadClass.metro]: COLORS.roads_metro.deck,
  [RoadClass.track]: COLORS.roads_unknown.deck,
  [RoadClass.other]: COLORS.roads_unknown.deck,
};
const roadLineSize: Record<RoadClass, ScaleLevel> = {
  [RoadClass.class_a]: 0,
  [RoadClass.class_b]: 1,
  [RoadClass.class_c]: 2,
  [RoadClass.metro]: 2,
  [RoadClass.track]: 2,
  [RoadClass.other]: 2,
};
function roadsViewLayer(asset_id) {
  const roadClass = roadClassLookup[asset_id];
  return infrastructureViewLayer(asset_id, ({ zoom, styleParams }) => [
    { minZoom: NETWORKS_METADATA[asset_id].minZoom ?? 4 },
    strokeColor(infraStyle(asset_id, roadColor[roadClass], styleParams)),
    lineStyle(zoom, roadLineSize[roadClass]),
  ]);
}

function potableNodesViewLayer(asset_id) {
  return infrastructureViewLayer(asset_id, ({ zoom, styleParams }) => [
    iconType('inv-triangle'),
    iconSize(zoom, 1),
    iconColor(infraStyle(asset_id, COLORS.water_supply.deck, styleParams)),
  ]);
}

function wastewaterNodesViewLayer(asset_id) {
  return infrastructureViewLayer(asset_id, ({ zoom, styleParams }) => [
    iconType('inv-triangle'),
    iconSize(zoom, 1),
    iconColor(infraStyle(asset_id, COLORS.water_wastewater.deck, styleParams)),
  ]);
}

export const INFRASTRUCTURE_VIEW_LAYERS = makeConfig<ViewLayer, string>([
  infrastructureViewLayer('elec_edges_high', ({ zoom, styleParams }) => [
    strokeColor(infraStyle('elec_edges_high', COLORS.electricity_high.deck, styleParams)),
    lineStyle(zoom, 1),
  ]),
  infrastructureViewLayer('elec_edges_low', ({ zoom, styleParams }) => [
    { minZoom: NETWORKS_METADATA.elec_edges_low.minZoom },
    strokeColor(infraStyle('elec_edges_low', COLORS.electricity_low.deck, styleParams)),
    lineStyle(zoom, 2),
  ]),

  infrastructureViewLayer('elec_nodes_diesel', ({ zoom, styleParams }) => [
    iconType('square'),
    iconColor(infraStyle('elec_nodes_diesel', COLORS.elec_nodes_diesel.deck, styleParams)),
    iconSize(zoom, 1),
  ]),
  infrastructureViewLayer('elec_nodes_gas', ({ zoom, styleParams }) => [
    iconType('square'),
    iconColor(infraStyle('elec_nodes_gas', COLORS.elec_nodes_gas.deck, styleParams)),
    iconSize(zoom, 1),
  ]),
  infrastructureViewLayer('elec_nodes_hydro', ({ zoom, styleParams }) => [
    iconType('square'),
    iconColor(infraStyle('elec_nodes_hydro', COLORS.elec_nodes_hydro.deck, styleParams)),
    iconSize(zoom, 1),
  ]),
  infrastructureViewLayer('elec_nodes_solar', ({ zoom, styleParams }) => [
    iconType('square'),
    iconColor(infraStyle('elec_nodes_solar', COLORS.elec_nodes_solar.deck, styleParams)),
    iconSize(zoom, 1),
  ]),
  infrastructureViewLayer('elec_nodes_wind', ({ zoom, styleParams }) => [
    iconType('square'),
    iconColor(infraStyle('elec_nodes_wind', COLORS.elec_nodes_wind.deck, styleParams)),
    iconSize(zoom, 1),
  ]),

  infrastructureViewLayer('elec_nodes_demand', ({ zoom, styleParams }) => [
    { minZoom: NETWORKS_METADATA.elec_nodes_demand.minZoom },
    fillColor(infraStyle('elec_nodes_demand', COLORS.electricity_demand.deck, styleParams)),
    pointRadius(zoom, 3),
    border(),
  ]),
  infrastructureViewLayer('elec_nodes_pole', ({ zoom, styleParams }) => [
    { minZoom: NETWORKS_METADATA.elec_nodes_pole.minZoom },
    fillColor(infraStyle('elec_nodes_pole', COLORS.electricity_low.deck, styleParams)),
    pointRadius(zoom, 3),
    border(),
  ]),
  infrastructureViewLayer('elec_nodes_substation', ({ zoom, styleParams }) => [
    border(),
    fillColor(infraStyle('elec_nodes_substation', COLORS.electricity_high.deck, styleParams)),
    pointRadius(zoom, 1),
  ]),

  infrastructureViewLayer('rail_edges', ({ zoom, styleParams }) => [
    strokeColor(infraStyle('rail_edges', COLORS.railway.deck, styleParams)),
    lineStyle(zoom, 1),
  ]),
  infrastructureViewLayer('rail_stations', ({ zoom, styleParams }) => [
    border(),
    fillColor(infraStyle('rail_stations', COLORS.railway.deck, styleParams)),
    pointRadius(zoom, 1),
  ]),
  infrastructureViewLayer('rail_junctions', ({ zoom, styleParams }) => [
    iconType('diamond'),
    iconColor(infraStyle('rail_junctions', COLORS.railway.deck, styleParams)),
    iconSize(zoom, 2),
  ]),

  roadsViewLayer('road_edges_class_a'),
  roadsViewLayer('road_edges_class_b'),
  roadsViewLayer('road_edges_class_c'),
  roadsViewLayer('road_edges_metro'),
  roadsViewLayer('road_edges_track'),
  roadsViewLayer('road_edges_other'),

  infrastructureViewLayer('road_bridges', ({ zoom, styleParams }) => [
    iconType('diamond'),
    iconColor(infraStyle('road_bridges', COLORS.bridges.deck, styleParams)),
    iconSize(zoom, 1),
  ]),
  infrastructureViewLayer('airport_runways', ({ zoom, styleParams }) => [
    zoom >= 10 && border(),
    fillColor(infraStyle('airport_runways', COLORS.airport_runways.deck, styleParams)),
  ]),
  infrastructureViewLayer('airport_terminals', ({ zoom, styleParams }) => [
    zoom >= 10 && border(),
    fillColor(infraStyle('airport_terminals', COLORS.airport_terminals.deck, styleParams)),
  ]),
  infrastructureViewLayer('port_areas_break', ({ zoom, styleParams }) => [
    zoom >= 10 && border(),
    fillColor(infraStyle('port_areas_break', COLORS.port_areas_break.deck, styleParams)),
  ]),
  infrastructureViewLayer('port_areas_container', ({ zoom, styleParams }) => [
    zoom >= 10 && border(),
    fillColor(infraStyle('port_areas_container', COLORS.port_areas_container.deck, styleParams)),
  ]),
  infrastructureViewLayer('port_areas_industry', ({ zoom, styleParams }) => [
    zoom >= 10 && border(),
    fillColor(infraStyle('port_areas_industry', COLORS.port_areas_industry.deck, styleParams)),
  ]),
  infrastructureViewLayer('port_areas_silo', ({ zoom, styleParams }) => [
    zoom >= 10 && border(),
    fillColor(infraStyle('port_areas_silo', COLORS.port_areas_silo.deck, styleParams)),
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
    iconType('inv-triangle'),
    iconSize(zoom, 1),
    iconColor(infraStyle('water_irrigation_nodes', COLORS.water_irrigation.deck, styleParams)),
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
