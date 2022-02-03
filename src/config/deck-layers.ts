import { makeConfig } from 'lib/helpers';
import { border, lineStyle, pointRadius, vectorColor } from 'lib/deck-layers/utils';

import { COLORS } from './colors';
import { hazardDeckLayer } from './deck-layers/hazard-layer';
import { infrastructureLayer } from './deck-layers/infrastructure-layer';

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

export const DECK_LAYERS = makeConfig<any, string>([
  {
    id: 'elec_edges_high',
    spatialType: 'vector',
    fn: ({ props, zoom, styleParams, selectedFeatureId }) =>
      infrastructureLayer(
        { selectedFeatureId },
        props,
        {
          data: '/vector/data/elec_edges_high.json',
        },
        vectorColor('stroke', COLORS.electricity_high.deck, styleParams),
        lineStyle(zoom),
      ),
  },
  {
    id: 'elec_edges_low',
    spatialType: 'vector',
    fn: ({ props, zoom, styleParams, selectedFeatureId }) =>
      infrastructureLayer(
        { selectedFeatureId },
        props,
        {
          data: '/vector/data/elec_edges_low.json',
        },
        vectorColor('stroke', COLORS.electricity_low.deck, styleParams),
        lineStyle(zoom),
      ),
  },
  {
    id: 'elec_nodes_source',
    spatialType: 'vector',
    fn: ({ props, zoom, styleParams, selectedFeatureId }) =>
      infrastructureLayer(
        { selectedFeatureId },
        props,
        {
          data: '/vector/data/elec_nodes_source.json',
        },
        border(),
        vectorColor('fill', COLORS.electricity_high.deck, styleParams),
        pointRadius(zoom),
      ),
  },
  {
    id: 'elec_nodes_sink',
    spatialType: 'vector',
    fn: ({ props, zoom, styleParams, selectedFeatureId }) =>
      infrastructureLayer(
        { selectedFeatureId },
        props,
        {
          data: '/vector/data/elec_nodes_sink.json',
        },
        border(),
        vectorColor('fill', COLORS.electricity_low.deck, styleParams),
        pointRadius(zoom),
      ),
  },
  {
    id: 'elec_nodes_junction',
    spatialType: 'vector',
    fn: ({ props, zoom, styleParams, selectedFeatureId }) =>
      infrastructureLayer(
        { selectedFeatureId },
        props,
        {
          data: '/vector/data/elec_nodes_junction.json',
        },
        border(),
        vectorColor('fill', COLORS.electricity_unknown.deck, styleParams),
        pointRadius(zoom),
      ),
  },
  {
    id: 'rail_edges',
    spatialType: 'vector',
    fn: ({ props, zoom, styleParams, selectedFeatureId }) =>
      infrastructureLayer(
        { selectedFeatureId },
        props,
        {
          data: '/vector/data/rail_edges.json',
        },
        vectorColor('stroke', COLORS.railway.deck, styleParams),
        lineStyle(zoom),
      ),
  },
  {
    id: 'rail_nodes',
    spatialType: 'vector',
    fn: ({ props, zoom, styleParams, selectedFeatureId }) =>
      infrastructureLayer(
        { selectedFeatureId },
        props,
        {
          data: '/vector/data/rail_nodes.json',
        },
        border(),
        vectorColor('fill', COLORS.railway.deck, styleParams),
        pointRadius(zoom),
      ),
  },
  {
    id: 'road_edges',
    spatialType: 'vector',
    fn: ({ props, zoom, styleParams, selectedFeatureId }) =>
      infrastructureLayer(
        { selectedFeatureId },
        props,
        {
          data: '/vector/data/road_edges.json',
        },
        vectorColor('stroke', (x) => roadColor[roadClassLookup[x.properties.road_class]], styleParams),
        lineStyle(zoom),
      ),
  },
  {
    id: 'road_bridges',
    spatialType: 'vector',
    fn: ({ props, zoom, styleParams, selectedFeatureId }) =>
      infrastructureLayer(
        { selectedFeatureId },
        props,
        {
          data: '/vector/data/road_bridges.json',
        },
        border(),
        vectorColor('fill', COLORS.bridges.deck, styleParams),
        pointRadius(zoom),
      ),
  },
  {
    id: 'airport_areas',
    spatialType: 'vector',
    fn: ({ props, zoom, styleParams, selectedFeatureId }) =>
      infrastructureLayer(
        { selectedFeatureId },
        props,
        {
          data: '/vector/data/airport_areas.json',
        },
        border([0, 0, 0]),
        vectorColor('fill', COLORS.airports.deck, styleParams),
      ),
  },
  {
    id: 'port_areas',
    spatialType: 'vector',
    fn: ({ props, zoom, styleParams, selectedFeatureId }) =>
      infrastructureLayer(
        { selectedFeatureId },
        props,
        {
          data: '/vector/data/port_areas.json',
        },
        border(),
        vectorColor('fill', COLORS.ports.deck, styleParams),
      ),
  },
  {
    id: 'water_potable_edges',
    spatialType: 'vector',
    fn: ({ props, zoom, styleParams, selectedFeatureId }) =>
      infrastructureLayer(
        { selectedFeatureId },
        props,
        {
          data: '/vector/data/water_potable_edges.json',
        },
        lineStyle(zoom),
        vectorColor('stroke', COLORS.water_supply.deck, styleParams),
      ),
  },
  {
    id: 'water_potable_nodes',
    spatialType: 'vector',
    fn: ({ props, zoom, styleParams, selectedFeatureId }) =>
      infrastructureLayer(
        { selectedFeatureId },
        props,
        {
          data: '/vector/data/water_potable_nodes.json',
        },
        border(),
        pointRadius(zoom),
        vectorColor('fill', COLORS.water_supply.deck, styleParams),
      ),
  },
  {
    id: 'water_irrigation_edges',
    spatialType: 'vector',
    fn: ({ props, zoom, styleParams, selectedFeatureId }) =>
      infrastructureLayer(
        { selectedFeatureId },
        props,
        {
          data: '/vector/data/water_irrigation_edges.json',
        },
        lineStyle(zoom),
        vectorColor('stroke', COLORS.water_irrigation.deck, styleParams),
      ),
  },
  {
    id: 'water_irrigation_nodes',
    spatialType: 'vector',
    fn: ({ props, zoom, styleParams, selectedFeatureId }) =>
      infrastructureLayer(
        { selectedFeatureId },
        props,
        {
          data: '/vector/data/water_irrigation_nodes.json',
        },
        border(),
        pointRadius(zoom),
        vectorColor('fill', COLORS.water_irrigation.deck, styleParams),
      ),
  },
  {
    id: 'water_waste_edges',
    spatialType: 'vector',
    fn: ({ props, zoom, styleParams, selectedFeatureId }) =>
      infrastructureLayer(
        { selectedFeatureId },
        props,
        {
          data: '/vector/data/water_waste_edges.json',
        },
        lineStyle(zoom),

        vectorColor('stroke', COLORS.water_wastewater.deck, styleParams),
      ),
  },
  {
    id: 'water_waste_nodes',
    spatialType: 'vector',
    fn: ({ props, zoom, styleParams, selectedFeatureId }) =>
      infrastructureLayer(
        { selectedFeatureId },
        props,
        {
          data: '/vector/data/water_waste_nodes.json',
        },
        border(),
        pointRadius(zoom),
        vectorColor('fill', COLORS.water_wastewater.deck, styleParams),
      ),
  },

  hazardDeckLayer('fluvial', 20, 'baseline', 2010, 'None'),
  hazardDeckLayer('fluvial', 50, 'baseline', 2010, 'None'),
  hazardDeckLayer('fluvial', 100, 'baseline', 2010, 'None'),
  hazardDeckLayer('fluvial', 200, 'baseline', 2010, 'None'),
  hazardDeckLayer('fluvial', 500, 'baseline', 2010, 'None'),
  hazardDeckLayer('fluvial', 1500, 'baseline', 2010, 'None'),

  hazardDeckLayer('surface', 20, 'baseline', 2010, 'None'),
  hazardDeckLayer('surface', 50, 'baseline', 2010, 'None'),
  hazardDeckLayer('surface', 100, 'baseline', 2010, 'None'),
  hazardDeckLayer('surface', 200, 'baseline', 2010, 'None'),
  hazardDeckLayer('surface', 500, 'baseline', 2010, 'None'),
  hazardDeckLayer('surface', 1500, 'baseline', 2010, 'None'),

  hazardDeckLayer('coastal', 1, '2.6', 2050, 'None'),
  hazardDeckLayer('coastal', 2, '2.6', 2050, 'None'),
  hazardDeckLayer('coastal', 5, '2.6', 2050, 'None'),
  hazardDeckLayer('coastal', 10, '2.6', 2050, 'None'),
  hazardDeckLayer('coastal', 50, '2.6', 2050, 'None'),
  hazardDeckLayer('coastal', 100, '2.6', 2050, 'None'),
  hazardDeckLayer('coastal', 1, '2.6', 2100, 'None'),
  hazardDeckLayer('coastal', 2, '2.6', 2100, 'None'),
  hazardDeckLayer('coastal', 5, '2.6', 2100, 'None'),
  hazardDeckLayer('coastal', 10, '2.6', 2100, 'None'),
  hazardDeckLayer('coastal', 50, '2.6', 2100, 'None'),
  hazardDeckLayer('coastal', 100, '2.6', 2100, 'None'),
  hazardDeckLayer('coastal', 1, 'baseline', 2010, 'None'),
  hazardDeckLayer('coastal', 2, 'baseline', 2010, 'None'),
  hazardDeckLayer('coastal', 5, 'baseline', 2010, 'None'),
  hazardDeckLayer('coastal', 10, 'baseline', 2010, 'None'),
  hazardDeckLayer('coastal', 50, 'baseline', 2010, 'None'),
  hazardDeckLayer('coastal', 100, 'baseline', 2010, 'None'),
  hazardDeckLayer('coastal', 1, '4.5', 2030, 'None'),
  hazardDeckLayer('coastal', 2, '4.5', 2030, 'None'),
  hazardDeckLayer('coastal', 5, '4.5', 2030, 'None'),
  hazardDeckLayer('coastal', 10, '4.5', 2030, 'None'),
  hazardDeckLayer('coastal', 50, '4.5', 2030, 'None'),
  hazardDeckLayer('coastal', 100, '4.5', 2030, 'None'),
  hazardDeckLayer('coastal', 1, '4.5', 2050, 'None'),
  hazardDeckLayer('coastal', 2, '4.5', 2050, 'None'),
  hazardDeckLayer('coastal', 5, '4.5', 2050, 'None'),
  hazardDeckLayer('coastal', 10, '4.5', 2050, 'None'),
  hazardDeckLayer('coastal', 50, '4.5', 2050, 'None'),
  hazardDeckLayer('coastal', 100, '4.5', 2050, 'None'),
  hazardDeckLayer('coastal', 1, '4.5', 2070, 'None'),
  hazardDeckLayer('coastal', 2, '4.5', 2070, 'None'),
  hazardDeckLayer('coastal', 5, '4.5', 2070, 'None'),
  hazardDeckLayer('coastal', 10, '4.5', 2070, 'None'),
  hazardDeckLayer('coastal', 50, '4.5', 2070, 'None'),
  hazardDeckLayer('coastal', 100, '4.5', 2070, 'None'),
  hazardDeckLayer('coastal', 1, '4.5', 2100, 'None'),
  hazardDeckLayer('coastal', 2, '4.5', 2100, 'None'),
  hazardDeckLayer('coastal', 5, '4.5', 2100, 'None'),
  hazardDeckLayer('coastal', 10, '4.5', 2100, 'None'),
  hazardDeckLayer('coastal', 50, '4.5', 2100, 'None'),
  hazardDeckLayer('coastal', 100, '4.5', 2100, 'None'),
  hazardDeckLayer('coastal', 1, '8.5', 2030, 'None'),
  hazardDeckLayer('coastal', 2, '8.5', 2030, 'None'),
  hazardDeckLayer('coastal', 5, '8.5', 2030, 'None'),
  hazardDeckLayer('coastal', 10, '8.5', 2030, 'None'),
  hazardDeckLayer('coastal', 50, '8.5', 2030, 'None'),
  hazardDeckLayer('coastal', 100, '8.5', 2030, 'None'),
  hazardDeckLayer('coastal', 1, '8.5', 2050, 'None'),
  hazardDeckLayer('coastal', 2, '8.5', 2050, 'None'),
  hazardDeckLayer('coastal', 5, '8.5', 2050, 'None'),
  hazardDeckLayer('coastal', 10, '8.5', 2050, 'None'),
  hazardDeckLayer('coastal', 50, '8.5', 2050, 'None'),
  hazardDeckLayer('coastal', 100, '8.5', 2050, 'None'),
  hazardDeckLayer('coastal', 1, '8.5', 2070, 'None'),
  hazardDeckLayer('coastal', 2, '8.5', 2070, 'None'),
  hazardDeckLayer('coastal', 5, '8.5', 2070, 'None'),
  hazardDeckLayer('coastal', 10, '8.5', 2070, 'None'),
  hazardDeckLayer('coastal', 50, '8.5', 2070, 'None'),
  hazardDeckLayer('coastal', 100, '8.5', 2070, 'None'),
  hazardDeckLayer('coastal', 1, '8.5', 2100, 'None'),
  hazardDeckLayer('coastal', 2, '8.5', 2100, 'None'),
  hazardDeckLayer('coastal', 5, '8.5', 2100, 'None'),
  hazardDeckLayer('coastal', 10, '8.5', 2100, 'None'),
  hazardDeckLayer('coastal', 50, '8.5', 2100, 'None'),
  hazardDeckLayer('coastal', 100, '8.5', 2100, 'None'),

  hazardDeckLayer('cyclone', 10, '4.5', 2050, 5),
  hazardDeckLayer('cyclone', 10, '4.5', 2050, 50),
  hazardDeckLayer('cyclone', 10, '4.5', 2050, 95),
  hazardDeckLayer('cyclone', 20, '4.5', 2050, 5),
  hazardDeckLayer('cyclone', 20, '4.5', 2050, 50),
  hazardDeckLayer('cyclone', 20, '4.5', 2050, 95),
  hazardDeckLayer('cyclone', 30, '4.5', 2050, 5),
  hazardDeckLayer('cyclone', 30, '4.5', 2050, 50),
  hazardDeckLayer('cyclone', 30, '4.5', 2050, 95),
  hazardDeckLayer('cyclone', 40, '4.5', 2050, 5),
  hazardDeckLayer('cyclone', 40, '4.5', 2050, 50),
  hazardDeckLayer('cyclone', 40, '4.5', 2050, 95),
  hazardDeckLayer('cyclone', 50, '4.5', 2050, 5),
  hazardDeckLayer('cyclone', 50, '4.5', 2050, 50),
  hazardDeckLayer('cyclone', 50, '4.5', 2050, 95),
  hazardDeckLayer('cyclone', 60, '4.5', 2050, 5),
  hazardDeckLayer('cyclone', 60, '4.5', 2050, 50),
  hazardDeckLayer('cyclone', 60, '4.5', 2050, 95),
  hazardDeckLayer('cyclone', 70, '4.5', 2050, 5),
  hazardDeckLayer('cyclone', 70, '4.5', 2050, 50),
  hazardDeckLayer('cyclone', 70, '4.5', 2050, 95),
  hazardDeckLayer('cyclone', 80, '4.5', 2050, 5),
  hazardDeckLayer('cyclone', 80, '4.5', 2050, 50),
  hazardDeckLayer('cyclone', 80, '4.5', 2050, 95),
  hazardDeckLayer('cyclone', 90, '4.5', 2050, 5),
  hazardDeckLayer('cyclone', 90, '4.5', 2050, 50),
  hazardDeckLayer('cyclone', 90, '4.5', 2050, 95),
  hazardDeckLayer('cyclone', 100, '4.5', 2050, 5),
  hazardDeckLayer('cyclone', 100, '4.5', 2050, 50),
  hazardDeckLayer('cyclone', 100, '4.5', 2050, 95),
  hazardDeckLayer('cyclone', 200, '4.5', 2050, 5),
  hazardDeckLayer('cyclone', 200, '4.5', 2050, 50),
  hazardDeckLayer('cyclone', 200, '4.5', 2050, 95),
  hazardDeckLayer('cyclone', 300, '4.5', 2050, 5),
  hazardDeckLayer('cyclone', 300, '4.5', 2050, 50),
  hazardDeckLayer('cyclone', 300, '4.5', 2050, 95),
  hazardDeckLayer('cyclone', 400, '4.5', 2050, 5),
  hazardDeckLayer('cyclone', 400, '4.5', 2050, 50),
  hazardDeckLayer('cyclone', 400, '4.5', 2050, 95),
  hazardDeckLayer('cyclone', 500, '4.5', 2050, 5),
  hazardDeckLayer('cyclone', 500, '4.5', 2050, 50),
  hazardDeckLayer('cyclone', 500, '4.5', 2050, 95),
  hazardDeckLayer('cyclone', 600, '4.5', 2050, 5),
  hazardDeckLayer('cyclone', 600, '4.5', 2050, 50),
  hazardDeckLayer('cyclone', 600, '4.5', 2050, 95),
  hazardDeckLayer('cyclone', 700, '4.5', 2050, 5),
  hazardDeckLayer('cyclone', 700, '4.5', 2050, 50),
  hazardDeckLayer('cyclone', 700, '4.5', 2050, 95),
  hazardDeckLayer('cyclone', 800, '4.5', 2050, 5),
  hazardDeckLayer('cyclone', 800, '4.5', 2050, 50),
  hazardDeckLayer('cyclone', 800, '4.5', 2050, 95),
  hazardDeckLayer('cyclone', 900, '4.5', 2050, 5),
  hazardDeckLayer('cyclone', 900, '4.5', 2050, 50),
  hazardDeckLayer('cyclone', 900, '4.5', 2050, 95),
  hazardDeckLayer('cyclone', 1000, '4.5', 2050, 5),
  hazardDeckLayer('cyclone', 1000, '4.5', 2050, 50),
  hazardDeckLayer('cyclone', 1000, '4.5', 2050, 95),
  hazardDeckLayer('cyclone', 2000, '4.5', 2050, 5),
  hazardDeckLayer('cyclone', 2000, '4.5', 2050, 50),
  hazardDeckLayer('cyclone', 2000, '4.5', 2050, 95),
  hazardDeckLayer('cyclone', 3000, '4.5', 2050, 5),
  hazardDeckLayer('cyclone', 3000, '4.5', 2050, 50),
  hazardDeckLayer('cyclone', 3000, '4.5', 2050, 95),
  hazardDeckLayer('cyclone', 4000, '4.5', 2050, 5),
  hazardDeckLayer('cyclone', 4000, '4.5', 2050, 50),
  hazardDeckLayer('cyclone', 4000, '4.5', 2050, 95),
  hazardDeckLayer('cyclone', 5000, '4.5', 2050, 5),
  hazardDeckLayer('cyclone', 5000, '4.5', 2050, 50),
  hazardDeckLayer('cyclone', 5000, '4.5', 2050, 95),
  hazardDeckLayer('cyclone', 6000, '4.5', 2050, 5),
  hazardDeckLayer('cyclone', 6000, '4.5', 2050, 50),
  hazardDeckLayer('cyclone', 6000, '4.5', 2050, 95),
  hazardDeckLayer('cyclone', 7000, '4.5', 2050, 5),
  hazardDeckLayer('cyclone', 7000, '4.5', 2050, 50),
  hazardDeckLayer('cyclone', 7000, '4.5', 2050, 95),
  hazardDeckLayer('cyclone', 8000, '4.5', 2050, 5),
  hazardDeckLayer('cyclone', 8000, '4.5', 2050, 50),
  hazardDeckLayer('cyclone', 8000, '4.5', 2050, 95),
  hazardDeckLayer('cyclone', 9000, '4.5', 2050, 5),
  hazardDeckLayer('cyclone', 9000, '4.5', 2050, 50),
  hazardDeckLayer('cyclone', 9000, '4.5', 2050, 95),
  hazardDeckLayer('cyclone', 10000, '4.5', 2050, 5),
  hazardDeckLayer('cyclone', 10000, '4.5', 2050, 50),
  hazardDeckLayer('cyclone', 10000, '4.5', 2050, 95),
  hazardDeckLayer('cyclone', 10, '4.5', 2100, 5),
  hazardDeckLayer('cyclone', 10, '4.5', 2100, 50),
  hazardDeckLayer('cyclone', 10, '4.5', 2100, 95),
  hazardDeckLayer('cyclone', 20, '4.5', 2100, 5),
  hazardDeckLayer('cyclone', 20, '4.5', 2100, 50),
  hazardDeckLayer('cyclone', 20, '4.5', 2100, 95),
  hazardDeckLayer('cyclone', 30, '4.5', 2100, 5),
  hazardDeckLayer('cyclone', 30, '4.5', 2100, 50),
  hazardDeckLayer('cyclone', 30, '4.5', 2100, 95),
  hazardDeckLayer('cyclone', 40, '4.5', 2100, 5),
  hazardDeckLayer('cyclone', 40, '4.5', 2100, 50),
  hazardDeckLayer('cyclone', 40, '4.5', 2100, 95),
  hazardDeckLayer('cyclone', 50, '4.5', 2100, 5),
  hazardDeckLayer('cyclone', 50, '4.5', 2100, 50),
  hazardDeckLayer('cyclone', 50, '4.5', 2100, 95),
  hazardDeckLayer('cyclone', 60, '4.5', 2100, 5),
  hazardDeckLayer('cyclone', 60, '4.5', 2100, 50),
  hazardDeckLayer('cyclone', 60, '4.5', 2100, 95),
  hazardDeckLayer('cyclone', 70, '4.5', 2100, 5),
  hazardDeckLayer('cyclone', 70, '4.5', 2100, 50),
  hazardDeckLayer('cyclone', 70, '4.5', 2100, 95),
  hazardDeckLayer('cyclone', 80, '4.5', 2100, 5),
  hazardDeckLayer('cyclone', 80, '4.5', 2100, 50),
  hazardDeckLayer('cyclone', 80, '4.5', 2100, 95),
  hazardDeckLayer('cyclone', 90, '4.5', 2100, 5),
  hazardDeckLayer('cyclone', 90, '4.5', 2100, 50),
  hazardDeckLayer('cyclone', 90, '4.5', 2100, 95),
  hazardDeckLayer('cyclone', 100, '4.5', 2100, 5),
  hazardDeckLayer('cyclone', 100, '4.5', 2100, 50),
  hazardDeckLayer('cyclone', 100, '4.5', 2100, 95),
  hazardDeckLayer('cyclone', 200, '4.5', 2100, 5),
  hazardDeckLayer('cyclone', 200, '4.5', 2100, 50),
  hazardDeckLayer('cyclone', 200, '4.5', 2100, 95),
  hazardDeckLayer('cyclone', 300, '4.5', 2100, 5),
  hazardDeckLayer('cyclone', 300, '4.5', 2100, 50),
  hazardDeckLayer('cyclone', 300, '4.5', 2100, 95),
  hazardDeckLayer('cyclone', 400, '4.5', 2100, 5),
  hazardDeckLayer('cyclone', 400, '4.5', 2100, 50),
  hazardDeckLayer('cyclone', 400, '4.5', 2100, 95),
  hazardDeckLayer('cyclone', 500, '4.5', 2100, 5),
  hazardDeckLayer('cyclone', 500, '4.5', 2100, 50),
  hazardDeckLayer('cyclone', 500, '4.5', 2100, 95),
  hazardDeckLayer('cyclone', 600, '4.5', 2100, 5),
  hazardDeckLayer('cyclone', 600, '4.5', 2100, 50),
  hazardDeckLayer('cyclone', 600, '4.5', 2100, 95),
  hazardDeckLayer('cyclone', 700, '4.5', 2100, 5),
  hazardDeckLayer('cyclone', 700, '4.5', 2100, 50),
  hazardDeckLayer('cyclone', 700, '4.5', 2100, 95),
  hazardDeckLayer('cyclone', 800, '4.5', 2100, 5),
  hazardDeckLayer('cyclone', 800, '4.5', 2100, 50),
  hazardDeckLayer('cyclone', 800, '4.5', 2100, 95),
  hazardDeckLayer('cyclone', 900, '4.5', 2100, 5),
  hazardDeckLayer('cyclone', 900, '4.5', 2100, 50),
  hazardDeckLayer('cyclone', 900, '4.5', 2100, 95),
  hazardDeckLayer('cyclone', 1000, '4.5', 2100, 5),
  hazardDeckLayer('cyclone', 1000, '4.5', 2100, 50),
  hazardDeckLayer('cyclone', 1000, '4.5', 2100, 95),
  hazardDeckLayer('cyclone', 2000, '4.5', 2100, 5),
  hazardDeckLayer('cyclone', 2000, '4.5', 2100, 50),
  hazardDeckLayer('cyclone', 2000, '4.5', 2100, 95),
  hazardDeckLayer('cyclone', 3000, '4.5', 2100, 5),
  hazardDeckLayer('cyclone', 3000, '4.5', 2100, 50),
  hazardDeckLayer('cyclone', 3000, '4.5', 2100, 95),
  hazardDeckLayer('cyclone', 4000, '4.5', 2100, 5),
  hazardDeckLayer('cyclone', 4000, '4.5', 2100, 50),
  hazardDeckLayer('cyclone', 4000, '4.5', 2100, 95),
  hazardDeckLayer('cyclone', 5000, '4.5', 2100, 5),
  hazardDeckLayer('cyclone', 5000, '4.5', 2100, 50),
  hazardDeckLayer('cyclone', 5000, '4.5', 2100, 95),
  hazardDeckLayer('cyclone', 6000, '4.5', 2100, 5),
  hazardDeckLayer('cyclone', 6000, '4.5', 2100, 50),
  hazardDeckLayer('cyclone', 6000, '4.5', 2100, 95),
  hazardDeckLayer('cyclone', 7000, '4.5', 2100, 5),
  hazardDeckLayer('cyclone', 7000, '4.5', 2100, 50),
  hazardDeckLayer('cyclone', 7000, '4.5', 2100, 95),
  hazardDeckLayer('cyclone', 8000, '4.5', 2100, 5),
  hazardDeckLayer('cyclone', 8000, '4.5', 2100, 50),
  hazardDeckLayer('cyclone', 8000, '4.5', 2100, 95),
  hazardDeckLayer('cyclone', 9000, '4.5', 2100, 5),
  hazardDeckLayer('cyclone', 9000, '4.5', 2100, 50),
  hazardDeckLayer('cyclone', 9000, '4.5', 2100, 95),
  hazardDeckLayer('cyclone', 10000, '4.5', 2100, 5),
  hazardDeckLayer('cyclone', 10000, '4.5', 2100, 50),
  hazardDeckLayer('cyclone', 10000, '4.5', 2100, 95),
  hazardDeckLayer('cyclone', 10, '8.5', 2050, 5),
  hazardDeckLayer('cyclone', 10, '8.5', 2050, 50),
  hazardDeckLayer('cyclone', 10, '8.5', 2050, 95),
  hazardDeckLayer('cyclone', 20, '8.5', 2050, 5),
  hazardDeckLayer('cyclone', 20, '8.5', 2050, 50),
  hazardDeckLayer('cyclone', 20, '8.5', 2050, 95),
  hazardDeckLayer('cyclone', 30, '8.5', 2050, 5),
  hazardDeckLayer('cyclone', 30, '8.5', 2050, 50),
  hazardDeckLayer('cyclone', 30, '8.5', 2050, 95),
  hazardDeckLayer('cyclone', 40, '8.5', 2050, 5),
  hazardDeckLayer('cyclone', 40, '8.5', 2050, 50),
  hazardDeckLayer('cyclone', 40, '8.5', 2050, 95),
  hazardDeckLayer('cyclone', 50, '8.5', 2050, 5),
  hazardDeckLayer('cyclone', 50, '8.5', 2050, 50),
  hazardDeckLayer('cyclone', 50, '8.5', 2050, 95),
  hazardDeckLayer('cyclone', 60, '8.5', 2050, 5),
  hazardDeckLayer('cyclone', 60, '8.5', 2050, 50),
  hazardDeckLayer('cyclone', 60, '8.5', 2050, 95),
  hazardDeckLayer('cyclone', 70, '8.5', 2050, 5),
  hazardDeckLayer('cyclone', 70, '8.5', 2050, 50),
  hazardDeckLayer('cyclone', 70, '8.5', 2050, 95),
  hazardDeckLayer('cyclone', 80, '8.5', 2050, 5),
  hazardDeckLayer('cyclone', 80, '8.5', 2050, 50),
  hazardDeckLayer('cyclone', 80, '8.5', 2050, 95),
  hazardDeckLayer('cyclone', 90, '8.5', 2050, 5),
  hazardDeckLayer('cyclone', 90, '8.5', 2050, 50),
  hazardDeckLayer('cyclone', 90, '8.5', 2050, 95),
  hazardDeckLayer('cyclone', 100, '8.5', 2050, 5),
  hazardDeckLayer('cyclone', 100, '8.5', 2050, 50),
  hazardDeckLayer('cyclone', 100, '8.5', 2050, 95),
  hazardDeckLayer('cyclone', 200, '8.5', 2050, 5),
  hazardDeckLayer('cyclone', 200, '8.5', 2050, 50),
  hazardDeckLayer('cyclone', 200, '8.5', 2050, 95),
  hazardDeckLayer('cyclone', 300, '8.5', 2050, 5),
  hazardDeckLayer('cyclone', 300, '8.5', 2050, 50),
  hazardDeckLayer('cyclone', 300, '8.5', 2050, 95),
  hazardDeckLayer('cyclone', 400, '8.5', 2050, 5),
  hazardDeckLayer('cyclone', 400, '8.5', 2050, 50),
  hazardDeckLayer('cyclone', 400, '8.5', 2050, 95),
  hazardDeckLayer('cyclone', 500, '8.5', 2050, 5),
  hazardDeckLayer('cyclone', 500, '8.5', 2050, 50),
  hazardDeckLayer('cyclone', 500, '8.5', 2050, 95),
  hazardDeckLayer('cyclone', 600, '8.5', 2050, 5),
  hazardDeckLayer('cyclone', 600, '8.5', 2050, 50),
  hazardDeckLayer('cyclone', 600, '8.5', 2050, 95),
  hazardDeckLayer('cyclone', 700, '8.5', 2050, 5),
  hazardDeckLayer('cyclone', 700, '8.5', 2050, 50),
  hazardDeckLayer('cyclone', 700, '8.5', 2050, 95),
  hazardDeckLayer('cyclone', 800, '8.5', 2050, 5),
  hazardDeckLayer('cyclone', 800, '8.5', 2050, 50),
  hazardDeckLayer('cyclone', 800, '8.5', 2050, 95),
  hazardDeckLayer('cyclone', 900, '8.5', 2050, 5),
  hazardDeckLayer('cyclone', 900, '8.5', 2050, 50),
  hazardDeckLayer('cyclone', 900, '8.5', 2050, 95),
  hazardDeckLayer('cyclone', 1000, '8.5', 2050, 5),
  hazardDeckLayer('cyclone', 1000, '8.5', 2050, 50),
  hazardDeckLayer('cyclone', 1000, '8.5', 2050, 95),
  hazardDeckLayer('cyclone', 2000, '8.5', 2050, 5),
  hazardDeckLayer('cyclone', 2000, '8.5', 2050, 50),
  hazardDeckLayer('cyclone', 2000, '8.5', 2050, 95),
  hazardDeckLayer('cyclone', 3000, '8.5', 2050, 5),
  hazardDeckLayer('cyclone', 3000, '8.5', 2050, 50),
  hazardDeckLayer('cyclone', 3000, '8.5', 2050, 95),
  hazardDeckLayer('cyclone', 4000, '8.5', 2050, 5),
  hazardDeckLayer('cyclone', 4000, '8.5', 2050, 50),
  hazardDeckLayer('cyclone', 4000, '8.5', 2050, 95),
  hazardDeckLayer('cyclone', 5000, '8.5', 2050, 5),
  hazardDeckLayer('cyclone', 5000, '8.5', 2050, 50),
  hazardDeckLayer('cyclone', 5000, '8.5', 2050, 95),
  hazardDeckLayer('cyclone', 6000, '8.5', 2050, 5),
  hazardDeckLayer('cyclone', 6000, '8.5', 2050, 50),
  hazardDeckLayer('cyclone', 6000, '8.5', 2050, 95),
  hazardDeckLayer('cyclone', 7000, '8.5', 2050, 5),
  hazardDeckLayer('cyclone', 7000, '8.5', 2050, 50),
  hazardDeckLayer('cyclone', 7000, '8.5', 2050, 95),
  hazardDeckLayer('cyclone', 8000, '8.5', 2050, 5),
  hazardDeckLayer('cyclone', 8000, '8.5', 2050, 50),
  hazardDeckLayer('cyclone', 8000, '8.5', 2050, 95),
  hazardDeckLayer('cyclone', 9000, '8.5', 2050, 5),
  hazardDeckLayer('cyclone', 9000, '8.5', 2050, 50),
  hazardDeckLayer('cyclone', 9000, '8.5', 2050, 95),
  hazardDeckLayer('cyclone', 10000, '8.5', 2050, 5),
  hazardDeckLayer('cyclone', 10000, '8.5', 2050, 50),
  hazardDeckLayer('cyclone', 10000, '8.5', 2050, 95),
  hazardDeckLayer('cyclone', 10, '8.5', 2100, 5),
  hazardDeckLayer('cyclone', 10, '8.5', 2100, 50),
  hazardDeckLayer('cyclone', 10, '8.5', 2100, 95),
  hazardDeckLayer('cyclone', 20, '8.5', 2100, 5),
  hazardDeckLayer('cyclone', 20, '8.5', 2100, 50),
  hazardDeckLayer('cyclone', 20, '8.5', 2100, 95),
  hazardDeckLayer('cyclone', 30, '8.5', 2100, 5),
  hazardDeckLayer('cyclone', 30, '8.5', 2100, 50),
  hazardDeckLayer('cyclone', 30, '8.5', 2100, 95),
  hazardDeckLayer('cyclone', 40, '8.5', 2100, 5),
  hazardDeckLayer('cyclone', 40, '8.5', 2100, 50),
  hazardDeckLayer('cyclone', 40, '8.5', 2100, 95),
  hazardDeckLayer('cyclone', 50, '8.5', 2100, 5),
  hazardDeckLayer('cyclone', 50, '8.5', 2100, 50),
  hazardDeckLayer('cyclone', 50, '8.5', 2100, 95),
  hazardDeckLayer('cyclone', 60, '8.5', 2100, 5),
  hazardDeckLayer('cyclone', 60, '8.5', 2100, 50),
  hazardDeckLayer('cyclone', 60, '8.5', 2100, 95),
  hazardDeckLayer('cyclone', 70, '8.5', 2100, 5),
  hazardDeckLayer('cyclone', 70, '8.5', 2100, 50),
  hazardDeckLayer('cyclone', 70, '8.5', 2100, 95),
  hazardDeckLayer('cyclone', 80, '8.5', 2100, 5),
  hazardDeckLayer('cyclone', 80, '8.5', 2100, 50),
  hazardDeckLayer('cyclone', 80, '8.5', 2100, 95),
  hazardDeckLayer('cyclone', 90, '8.5', 2100, 5),
  hazardDeckLayer('cyclone', 90, '8.5', 2100, 50),
  hazardDeckLayer('cyclone', 90, '8.5', 2100, 95),
  hazardDeckLayer('cyclone', 100, '8.5', 2100, 5),
  hazardDeckLayer('cyclone', 100, '8.5', 2100, 50),
  hazardDeckLayer('cyclone', 100, '8.5', 2100, 95),
  hazardDeckLayer('cyclone', 200, '8.5', 2100, 5),
  hazardDeckLayer('cyclone', 200, '8.5', 2100, 50),
  hazardDeckLayer('cyclone', 200, '8.5', 2100, 95),
  hazardDeckLayer('cyclone', 300, '8.5', 2100, 5),
  hazardDeckLayer('cyclone', 300, '8.5', 2100, 50),
  hazardDeckLayer('cyclone', 300, '8.5', 2100, 95),
  hazardDeckLayer('cyclone', 400, '8.5', 2100, 5),
  hazardDeckLayer('cyclone', 400, '8.5', 2100, 50),
  hazardDeckLayer('cyclone', 400, '8.5', 2100, 95),
  hazardDeckLayer('cyclone', 500, '8.5', 2100, 5),
  hazardDeckLayer('cyclone', 500, '8.5', 2100, 50),
  hazardDeckLayer('cyclone', 500, '8.5', 2100, 95),
  hazardDeckLayer('cyclone', 600, '8.5', 2100, 5),
  hazardDeckLayer('cyclone', 600, '8.5', 2100, 50),
  hazardDeckLayer('cyclone', 600, '8.5', 2100, 95),
  hazardDeckLayer('cyclone', 700, '8.5', 2100, 5),
  hazardDeckLayer('cyclone', 700, '8.5', 2100, 50),
  hazardDeckLayer('cyclone', 700, '8.5', 2100, 95),
  hazardDeckLayer('cyclone', 800, '8.5', 2100, 5),
  hazardDeckLayer('cyclone', 800, '8.5', 2100, 50),
  hazardDeckLayer('cyclone', 800, '8.5', 2100, 95),
  hazardDeckLayer('cyclone', 900, '8.5', 2100, 5),
  hazardDeckLayer('cyclone', 900, '8.5', 2100, 50),
  hazardDeckLayer('cyclone', 900, '8.5', 2100, 95),
  hazardDeckLayer('cyclone', 1000, '8.5', 2100, 5),
  hazardDeckLayer('cyclone', 1000, '8.5', 2100, 50),
  hazardDeckLayer('cyclone', 1000, '8.5', 2100, 95),
  hazardDeckLayer('cyclone', 2000, '8.5', 2100, 5),
  hazardDeckLayer('cyclone', 2000, '8.5', 2100, 50),
  hazardDeckLayer('cyclone', 2000, '8.5', 2100, 95),
  hazardDeckLayer('cyclone', 3000, '8.5', 2100, 5),
  hazardDeckLayer('cyclone', 3000, '8.5', 2100, 50),
  hazardDeckLayer('cyclone', 3000, '8.5', 2100, 95),
  hazardDeckLayer('cyclone', 4000, '8.5', 2100, 5),
  hazardDeckLayer('cyclone', 4000, '8.5', 2100, 50),
  hazardDeckLayer('cyclone', 4000, '8.5', 2100, 95),
  hazardDeckLayer('cyclone', 5000, '8.5', 2100, 5),
  hazardDeckLayer('cyclone', 5000, '8.5', 2100, 50),
  hazardDeckLayer('cyclone', 5000, '8.5', 2100, 95),
  hazardDeckLayer('cyclone', 6000, '8.5', 2100, 5),
  hazardDeckLayer('cyclone', 6000, '8.5', 2100, 50),
  hazardDeckLayer('cyclone', 6000, '8.5', 2100, 95),
  hazardDeckLayer('cyclone', 7000, '8.5', 2100, 5),
  hazardDeckLayer('cyclone', 7000, '8.5', 2100, 50),
  hazardDeckLayer('cyclone', 7000, '8.5', 2100, 95),
  hazardDeckLayer('cyclone', 8000, '8.5', 2100, 5),
  hazardDeckLayer('cyclone', 8000, '8.5', 2100, 50),
  hazardDeckLayer('cyclone', 8000, '8.5', 2100, 95),
  hazardDeckLayer('cyclone', 9000, '8.5', 2100, 5),
  hazardDeckLayer('cyclone', 9000, '8.5', 2100, 50),
  hazardDeckLayer('cyclone', 9000, '8.5', 2100, 95),
  hazardDeckLayer('cyclone', 10000, '8.5', 2100, 5),
  hazardDeckLayer('cyclone', 10000, '8.5', 2100, 50),
  hazardDeckLayer('cyclone', 10000, '8.5', 2100, 95),
  hazardDeckLayer('cyclone', 10, 'baseline', 2010, 5),
  hazardDeckLayer('cyclone', 10, 'baseline', 2010, 50),
  hazardDeckLayer('cyclone', 10, 'baseline', 2010, 95),
  hazardDeckLayer('cyclone', 20, 'baseline', 2010, 5),
  hazardDeckLayer('cyclone', 20, 'baseline', 2010, 50),
  hazardDeckLayer('cyclone', 20, 'baseline', 2010, 95),
  hazardDeckLayer('cyclone', 30, 'baseline', 2010, 5),
  hazardDeckLayer('cyclone', 30, 'baseline', 2010, 50),
  hazardDeckLayer('cyclone', 30, 'baseline', 2010, 95),
  hazardDeckLayer('cyclone', 40, 'baseline', 2010, 5),
  hazardDeckLayer('cyclone', 40, 'baseline', 2010, 50),
  hazardDeckLayer('cyclone', 40, 'baseline', 2010, 95),
  hazardDeckLayer('cyclone', 50, 'baseline', 2010, 5),
  hazardDeckLayer('cyclone', 50, 'baseline', 2010, 50),
  hazardDeckLayer('cyclone', 50, 'baseline', 2010, 95),
  hazardDeckLayer('cyclone', 60, 'baseline', 2010, 5),
  hazardDeckLayer('cyclone', 60, 'baseline', 2010, 50),
  hazardDeckLayer('cyclone', 60, 'baseline', 2010, 95),
  hazardDeckLayer('cyclone', 70, 'baseline', 2010, 5),
  hazardDeckLayer('cyclone', 70, 'baseline', 2010, 50),
  hazardDeckLayer('cyclone', 70, 'baseline', 2010, 95),
  hazardDeckLayer('cyclone', 80, 'baseline', 2010, 5),
  hazardDeckLayer('cyclone', 80, 'baseline', 2010, 50),
  hazardDeckLayer('cyclone', 80, 'baseline', 2010, 95),
  hazardDeckLayer('cyclone', 90, 'baseline', 2010, 5),
  hazardDeckLayer('cyclone', 90, 'baseline', 2010, 50),
  hazardDeckLayer('cyclone', 90, 'baseline', 2010, 95),
  hazardDeckLayer('cyclone', 100, 'baseline', 2010, 5),
  hazardDeckLayer('cyclone', 100, 'baseline', 2010, 50),
  hazardDeckLayer('cyclone', 100, 'baseline', 2010, 95),
  hazardDeckLayer('cyclone', 200, 'baseline', 2010, 5),
  hazardDeckLayer('cyclone', 200, 'baseline', 2010, 50),
  hazardDeckLayer('cyclone', 200, 'baseline', 2010, 95),
  hazardDeckLayer('cyclone', 300, 'baseline', 2010, 5),
  hazardDeckLayer('cyclone', 300, 'baseline', 2010, 50),
  hazardDeckLayer('cyclone', 300, 'baseline', 2010, 95),
  hazardDeckLayer('cyclone', 400, 'baseline', 2010, 5),
  hazardDeckLayer('cyclone', 400, 'baseline', 2010, 50),
  hazardDeckLayer('cyclone', 400, 'baseline', 2010, 95),
  hazardDeckLayer('cyclone', 500, 'baseline', 2010, 5),
  hazardDeckLayer('cyclone', 500, 'baseline', 2010, 50),
  hazardDeckLayer('cyclone', 500, 'baseline', 2010, 95),
  hazardDeckLayer('cyclone', 600, 'baseline', 2010, 5),
  hazardDeckLayer('cyclone', 600, 'baseline', 2010, 50),
  hazardDeckLayer('cyclone', 600, 'baseline', 2010, 95),
  hazardDeckLayer('cyclone', 700, 'baseline', 2010, 5),
  hazardDeckLayer('cyclone', 700, 'baseline', 2010, 50),
  hazardDeckLayer('cyclone', 700, 'baseline', 2010, 95),
  hazardDeckLayer('cyclone', 800, 'baseline', 2010, 5),
  hazardDeckLayer('cyclone', 800, 'baseline', 2010, 50),
  hazardDeckLayer('cyclone', 800, 'baseline', 2010, 95),
  hazardDeckLayer('cyclone', 900, 'baseline', 2010, 5),
  hazardDeckLayer('cyclone', 900, 'baseline', 2010, 50),
  hazardDeckLayer('cyclone', 900, 'baseline', 2010, 95),
  hazardDeckLayer('cyclone', 1000, 'baseline', 2010, 5),
  hazardDeckLayer('cyclone', 1000, 'baseline', 2010, 50),
  hazardDeckLayer('cyclone', 1000, 'baseline', 2010, 95),
  hazardDeckLayer('cyclone', 2000, 'baseline', 2010, 5),
  hazardDeckLayer('cyclone', 2000, 'baseline', 2010, 50),
  hazardDeckLayer('cyclone', 2000, 'baseline', 2010, 95),
  hazardDeckLayer('cyclone', 3000, 'baseline', 2010, 5),
  hazardDeckLayer('cyclone', 3000, 'baseline', 2010, 50),
  hazardDeckLayer('cyclone', 3000, 'baseline', 2010, 95),
  hazardDeckLayer('cyclone', 4000, 'baseline', 2010, 5),
  hazardDeckLayer('cyclone', 4000, 'baseline', 2010, 50),
  hazardDeckLayer('cyclone', 4000, 'baseline', 2010, 95),
  hazardDeckLayer('cyclone', 5000, 'baseline', 2010, 5),
  hazardDeckLayer('cyclone', 5000, 'baseline', 2010, 50),
  hazardDeckLayer('cyclone', 5000, 'baseline', 2010, 95),
  hazardDeckLayer('cyclone', 6000, 'baseline', 2010, 5),
  hazardDeckLayer('cyclone', 6000, 'baseline', 2010, 50),
  hazardDeckLayer('cyclone', 6000, 'baseline', 2010, 95),
  hazardDeckLayer('cyclone', 7000, 'baseline', 2010, 5),
  hazardDeckLayer('cyclone', 7000, 'baseline', 2010, 50),
  hazardDeckLayer('cyclone', 7000, 'baseline', 2010, 95),
  hazardDeckLayer('cyclone', 8000, 'baseline', 2010, 5),
  hazardDeckLayer('cyclone', 8000, 'baseline', 2010, 50),
  hazardDeckLayer('cyclone', 8000, 'baseline', 2010, 95),
  hazardDeckLayer('cyclone', 9000, 'baseline', 2010, 5),
  hazardDeckLayer('cyclone', 9000, 'baseline', 2010, 50),
  hazardDeckLayer('cyclone', 9000, 'baseline', 2010, 95),
  hazardDeckLayer('cyclone', 10000, 'baseline', 2010, 5),
  hazardDeckLayer('cyclone', 10000, 'baseline', 2010, 50),
  hazardDeckLayer('cyclone', 10000, 'baseline', 2010, 95),
]);
