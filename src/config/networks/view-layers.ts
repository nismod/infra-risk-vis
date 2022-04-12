import { makeConfig } from 'lib/helpers';
import { border, lineStyle, pointRadius } from 'lib/deck/props/style';

import { COLORS } from '../colors';
import { infrastructureViewLayer } from './infrastructure-view-layer';
import { ViewLayer } from 'lib/data-map/view-layers';
import { fillColor, strokeColor } from 'lib/deck/props/style';
import { dataColorMap } from 'lib/deck/props/color-map';
import { featureProperty } from 'lib/deck/props/data-source';
import { colorMapFromScheme } from 'config/color-maps';
import { getAssetDataAccessor } from 'config/assets/data-access';

enum RoadClass {
  class_a = 'class_a',
  class_b = 'class_b',
  class_c = 'class_c',
  metro = 'metro',
  other = 'other',
  track = 'track',
}

const roadClassLookup = {
  'CLASS A': RoadClass.class_a,
  'CLASS B': RoadClass.class_b,
  'CLASS C': RoadClass.class_c,
  METRO: RoadClass.metro,
  TRACK: RoadClass.track,
  OTHER: RoadClass.other,
};

const roadColor = {
  [RoadClass.class_a]: COLORS.roads_class_a.css,
  [RoadClass.class_b]: COLORS.roads_class_b.css,
  [RoadClass.class_c]: COLORS.roads_class_c.css,
  [RoadClass.metro]: COLORS.roads_class_metro.css,
  [RoadClass.track]: COLORS.roads_unknown.css,
  [RoadClass.other]: COLORS.roads_unknown.css,
};

function infraStyle(layer: string, defaultStyle, styleParams) {
  if (styleParams?.colorMap) {
    const { colorField, colorScheme } = styleParams.colorMap;
    return dataColorMap(getAssetDataAccessor(layer, colorField), colorMapFromScheme(colorScheme));
  }
  return defaultStyle;
}

const roadDefaultStyle = dataColorMap(featureProperty('road_class'), (x) => roadColor[roadClassLookup[x]]);

export const INFRASTRUCTURE_VIEW_LAYERS = makeConfig<ViewLayer, string>([
  infrastructureViewLayer('elec_edges_high', ({ zoom, styleParams }) => [
    strokeColor(infraStyle('elec_edges_high', COLORS.electricity_high.deck, styleParams)),
    lineStyle(zoom),
  ]),
  infrastructureViewLayer('elec_edges_low', ({ zoom, styleParams }) => [
    strokeColor(infraStyle('elec_edges_low', COLORS.electricity_low.deck, styleParams)),
    lineStyle(zoom),
  ]),
  infrastructureViewLayer('elec_nodes_source', ({ zoom, styleParams }) => [
    border(),
    fillColor(infraStyle('elec_nodes_source', COLORS.electricity_high.deck, styleParams)),
    pointRadius(zoom),
  ]),
  infrastructureViewLayer('elec_nodes_sink', ({ zoom, styleParams }) => [
    border(),
    fillColor(infraStyle('elec_nodes_sink', COLORS.electricity_low.deck, styleParams)),
    pointRadius(zoom),
  ]),
  infrastructureViewLayer('elec_nodes_junction', ({ zoom, styleParams }) => [
    border(),
    fillColor(infraStyle('elec_nodes_junction', COLORS.electricity_unknown.deck, styleParams)),
    pointRadius(zoom),
  ]),
  infrastructureViewLayer('rail_edges', ({ zoom, styleParams }) => [
    strokeColor(infraStyle('rail_edges', COLORS.railway.deck, styleParams)),
    lineStyle(zoom),
  ]),
  infrastructureViewLayer('rail_nodes', ({ zoom, styleParams }) => [
    border(),
    fillColor(infraStyle('rail_nodes', COLORS.railway.deck, styleParams)),
    pointRadius(zoom),
  ]),
  infrastructureViewLayer('road_edges', ({ zoom, styleParams }) => [
    strokeColor(infraStyle('road_edges', roadDefaultStyle, styleParams)),
    lineStyle(zoom),
  ]),
  infrastructureViewLayer('road_bridges', ({ zoom, styleParams }) => [
    border(),
    fillColor(infraStyle('road_bridges', COLORS.bridges.deck, styleParams)),
    pointRadius(zoom),
  ]),
  infrastructureViewLayer('airport_areas', ({ zoom, styleParams }) => [
    border([0, 0, 0]),
    fillColor(infraStyle('airport_areas', COLORS.airports.deck, styleParams)),
  ]),
  infrastructureViewLayer('port_areas', ({ zoom, styleParams }) => [
    border(),
    fillColor(infraStyle('port_areas', COLORS.ports.deck, styleParams)),
  ]),
  infrastructureViewLayer('water_potable_edges', ({ zoom, styleParams }) => [
    lineStyle(zoom),
    strokeColor(infraStyle('water_potable_edges', COLORS.water_supply.deck, styleParams)),
  ]),
  infrastructureViewLayer('water_potable_nodes', ({ zoom, styleParams }) => [
    border(),
    pointRadius(zoom),
    fillColor(infraStyle('water_potable_nodes', COLORS.water_supply.deck, styleParams)),
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
  infrastructureViewLayer('water_waste_edges', ({ zoom, styleParams }) => [
    lineStyle(zoom),
    strokeColor(infraStyle('water_waste_edges', COLORS.water_wastewater.deck, styleParams)),
  ]),
  infrastructureViewLayer('water_waste_nodes', ({ zoom, styleParams }) => [
    border(),
    pointRadius(zoom),
    fillColor(infraStyle('water_waste_nodes', COLORS.water_wastewater.deck, styleParams)),
  ]),
]);
