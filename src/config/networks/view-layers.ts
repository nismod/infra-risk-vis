import { makeConfig } from 'lib/helpers';
import { border, lineStyle, pointRadius } from 'lib/deck/props/style';

import { COLORS } from '../colors';
import { infrastructureViewLayer } from './infrastructure-view-layer';
import { ViewLayer } from 'lib/data-map/view-layers';
import { fillColor, strokeColor } from 'lib/deck/props/style';
import { dataColorMap } from 'lib/deck/props/color-map';
import { featureProperty } from 'lib/deck/props/data-source';
import { colorMapFromScheme } from 'config/color-maps';

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
  [RoadClass.class_a]: COLORS.roads_class_a.deck,
  [RoadClass.class_b]: COLORS.roads_class_b.deck,
  [RoadClass.class_c]: COLORS.roads_class_c.deck,
  [RoadClass.metro]: COLORS.roads_class_metro.deck,
  [RoadClass.track]: COLORS.roads_unknown.deck,
  [RoadClass.other]: COLORS.roads_unknown.deck,
};

function infraStyle(defaultStyle, styleParams) {
  if (styleParams?.colorMap) {
    const { colorField, colorScheme } = styleParams.colorMap;
    return dataColorMap(featureProperty(colorField), colorMapFromScheme(colorScheme));
  }
  return defaultStyle;
}

const roadDefaultStyle = dataColorMap(
  featureProperty('road_class'),
  (x) => roadColor[roadClassLookup[x]],
);

export const INFRASTRUCTURE_VIEW_LAYERS = makeConfig<ViewLayer, string>([
  infrastructureViewLayer('elec_edges_high', ({ zoom, styleParams }) => [
    strokeColor(infraStyle(COLORS.electricity_high.deck, styleParams)),
    lineStyle(zoom),
  ]),
  infrastructureViewLayer('elec_edges_low', ({ zoom, styleParams }) => [
    strokeColor(infraStyle(COLORS.electricity_low.deck, styleParams)),
    lineStyle(zoom),
  ]),
  infrastructureViewLayer('elec_nodes_source', ({ zoom, styleParams }) => [
    border(),
    fillColor(infraStyle(COLORS.electricity_high.deck, styleParams)),
    pointRadius(zoom),
  ]),
  infrastructureViewLayer('elec_nodes_sink', ({ zoom, styleParams }) => [
    border(),
    fillColor(infraStyle(COLORS.electricity_low.deck, styleParams)),
    pointRadius(zoom),
  ]),
  infrastructureViewLayer('elec_nodes_junction', ({ zoom, styleParams }) => [
    border(),
    fillColor(infraStyle(COLORS.electricity_unknown.deck, styleParams)),
    pointRadius(zoom),
  ]),
  infrastructureViewLayer('rail_edges', ({ zoom, styleParams }) => [
    strokeColor(infraStyle(COLORS.railway.deck, styleParams)),
    lineStyle(zoom),
  ]),
  infrastructureViewLayer('rail_nodes', ({ zoom, styleParams }) => [
    border(),
    fillColor(infraStyle(COLORS.railway.deck, styleParams)),
    pointRadius(zoom),
  ]),
  infrastructureViewLayer('road_edges', ({ zoom, styleParams }) => [
    strokeColor(infraStyle(roadDefaultStyle, styleParams)),
    lineStyle(zoom),
  ]),
  infrastructureViewLayer('road_bridges', ({ zoom, styleParams }) => [
    border(),
    fillColor(infraStyle(COLORS.bridges.deck, styleParams)),
    pointRadius(zoom),
  ]),
  infrastructureViewLayer('airport_areas', ({ zoom, styleParams }) => [
    border([0, 0, 0]),
    fillColor(infraStyle(COLORS.airports.deck, styleParams)),
  ]),
  infrastructureViewLayer('port_areas', ({ zoom, styleParams }) => [
    border(),
    fillColor(infraStyle(COLORS.ports.deck, styleParams)),
  ]),
  infrastructureViewLayer('water_potable_edges', ({ zoom, styleParams }) => [
    lineStyle(zoom),
    strokeColor(infraStyle(COLORS.water_supply.deck, styleParams)),
  ]),
  infrastructureViewLayer('water_potable_nodes', ({ zoom, styleParams }) => [
    border(),
    pointRadius(zoom),
    fillColor(infraStyle(COLORS.water_supply.deck, styleParams)),
  ]),
  infrastructureViewLayer('water_irrigation_edges', ({ zoom, styleParams }) => [
    lineStyle(zoom),
    strokeColor(infraStyle(COLORS.water_irrigation.deck, styleParams)),
  ]),
  infrastructureViewLayer('water_irrigation_nodes', ({ zoom, styleParams }) => [
    border(),
    pointRadius(zoom),
    fillColor(infraStyle(COLORS.water_irrigation.deck, styleParams)),
  ]),
  infrastructureViewLayer('water_waste_edges', ({ zoom, styleParams }) => [
    lineStyle(zoom),
    strokeColor(infraStyle(COLORS.water_wastewater.deck, styleParams)),
  ]),
  infrastructureViewLayer('water_waste_nodes', ({ zoom, styleParams }) => [
    border(),
    pointRadius(zoom),
    fillColor(infraStyle(COLORS.water_wastewater.deck, styleParams)),
  ]),
]);
