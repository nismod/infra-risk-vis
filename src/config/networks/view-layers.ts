import { makeConfig } from 'lib/helpers';
import { border, lineStyle, pointRadius, vectorColor } from 'lib/deck-layers/utils';

import { COLORS } from '../colors';
import { infrastructureViewLayer } from './infrastructure-view-layer';
import { ViewLayer } from 'lib/data-map/view-layers';

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

export const INFRASTRUCTURE_VIEW_LAYERS = makeConfig<ViewLayer, string>([
  infrastructureViewLayer('elec_edges_high', ({ zoom, styleParams }) => [
    vectorColor('stroke', COLORS.electricity_high.deck, styleParams),
    lineStyle(zoom),
  ]),
  infrastructureViewLayer('elec_edges_low', ({ zoom, styleParams }) => [
    vectorColor('stroke', COLORS.electricity_low.deck, styleParams),
    lineStyle(zoom),
  ]),
  infrastructureViewLayer('elec_nodes_source', ({ zoom, styleParams }) => [
    border(),
    vectorColor('fill', COLORS.electricity_high.deck, styleParams),
    pointRadius(zoom),
  ]),
  infrastructureViewLayer('elec_nodes_sink', ({ zoom, styleParams }) => [
    border(),
    vectorColor('fill', COLORS.electricity_low.deck, styleParams),
    pointRadius(zoom),
  ]),
  infrastructureViewLayer('elec_nodes_junction', ({ zoom, styleParams }) => [
    border(),
    vectorColor('fill', COLORS.electricity_unknown.deck, styleParams),
    pointRadius(zoom),
  ]),
  infrastructureViewLayer('rail_edges', ({ zoom, styleParams }) => [
    vectorColor('stroke', COLORS.railway.deck, styleParams),
    lineStyle(zoom),
  ]),
  infrastructureViewLayer('rail_nodes', ({ zoom, styleParams }) => [
    border(),
    vectorColor('fill', COLORS.railway.deck, styleParams),
    pointRadius(zoom),
  ]),
  infrastructureViewLayer('road_edges', ({ zoom, styleParams }) => [
    vectorColor('stroke', (x) => roadColor[roadClassLookup[x.properties.road_class]], styleParams),
    lineStyle(zoom),
  ]),
  infrastructureViewLayer('road_bridges', ({ zoom, styleParams }) => [
    border(),
    vectorColor('fill', COLORS.bridges.deck, styleParams),
    pointRadius(zoom),
  ]),
  infrastructureViewLayer('airport_areas', ({ zoom, styleParams }) => [
    border([0, 0, 0]),
    vectorColor('fill', COLORS.airports.deck, styleParams),
  ]),
  infrastructureViewLayer('port_areas', ({ zoom, styleParams }) => [
    border(),
    vectorColor('fill', COLORS.ports.deck, styleParams),
  ]),
  infrastructureViewLayer('water_potable_edges', ({ zoom, styleParams }) => [
    lineStyle(zoom),
    vectorColor('stroke', COLORS.water_supply.deck, styleParams),
  ]),
  infrastructureViewLayer('water_potable_nodes', ({ zoom, styleParams }) => [
    border(),
    pointRadius(zoom),
    vectorColor('fill', COLORS.water_supply.deck, styleParams),
  ]),
  infrastructureViewLayer('water_irrigation_edges', ({ zoom, styleParams }) => [
    lineStyle(zoom),
    vectorColor('stroke', COLORS.water_irrigation.deck, styleParams),
  ]),
  infrastructureViewLayer('water_irrigation_nodes', ({ zoom, styleParams }) => [
    border(),
    pointRadius(zoom),
    vectorColor('fill', COLORS.water_irrigation.deck, styleParams),
  ]),
  infrastructureViewLayer('water_waste_edges', ({ zoom, styleParams }) => [
    lineStyle(zoom),
    vectorColor('stroke', COLORS.water_wastewater.deck, styleParams),
  ]),
  infrastructureViewLayer('water_waste_nodes', ({ zoom, styleParams }) => [
    border(),
    pointRadius(zoom),
    vectorColor('fill', COLORS.water_wastewater.deck, styleParams),
  ]),
]);
