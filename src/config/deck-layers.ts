import { MVTLayer, TileLayer, BitmapLayer } from 'deck.gl';
import GL from '@luma.gl/constants';
import { DataFilterExtension } from '@deck.gl/extensions';

import { COLORS } from './colors';
import { makeConfig } from '../helpers';
import { getHazardId } from './layers';

function makeColormap(attr: string, valueColorMap: any) {
  return (x) => valueColorMap[x.properties[attr]];
}

const lineStyle = (zoom) => ({
  getLineWidth: 15,
  lineWidthUnit: 'meters',
  lineWidthMinPixels: 1,
  lineWidthMaxPixels: 5,
  lineJointRounded: true,
  lineCapRounded: true,

  // widthScale: 2 ** (15 - zoom),
});

const pointRadius = (zoom) => ({
  getPointRadius: 20,
  pointRadiusUnit: 'meters',
  pointRadiusMinPixels: 3,
  pointRadiusMaxPixels: 10,
  // radiusScale: 2 ** (15 - zoom),
});

function makeValueMapping<SourceValueType, ObjectType = any>(
  accessor: (x: ObjectType) => SourceValueType,
  sourceValues: SourceValueType[],
) {
  return <TargetValueType>(targetValues: TargetValueType[], defaultValue: TargetValueType) => {
    if (targetValues.length !== sourceValues.length) {
      throw new Error('Target domain length is not equal to the source domain length');
    }
    return (obj: ObjectType) => {
      const val = accessor(obj);
      const index = sourceValues.indexOf(val);
      if (index === -1) {
        return defaultValue;
      } else {
        return targetValues[index];
      }
    };
  };
}

const rasterColormaps = {
  fluvial: 'blues',
  coastal: 'greens',
  surface: 'purples',
};

const rasterColormapRanges = {
  fluvial: '[0,10]',
  coastal: '[0,3.5]',
  surface: '[0,10]',
};

function getBoundsForTile(tileProps) {
  const {
    bbox: { west, south, east, north },
  } = tileProps;

  return [west, south, east, north];
}

enum ElecVoltage {
  elec_edges_high = 'elec_edges_high',
  elec_edges_low = 'elec_edges_low',
}

const elecVoltageLookup = {
  'High Voltage': ElecVoltage.elec_edges_high,
  'Low Voltage': ElecVoltage.elec_edges_low,
};

const electricityColor = {
  [ElecVoltage.elec_edges_high]: COLORS.electricity_high.deck,
  [ElecVoltage.elec_edges_low]: COLORS.electricity_low.deck,
};

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
    id: 'elec_edges',
    type: 'MVTLayer',
    spatialType: 'vector',
    fn: ({ props, zoom, visibility }) =>
      new MVTLayer(props, {
        data: 'http://localhost:8080/data/elec_edges.json',
        dataTransform: (data) => {
          for (const objectProperty of data.lines.properties) {
            objectProperty.__logicalLayer = elecVoltageLookup[objectProperty.asset_type];
          }
          return data;
        },
        getLineColor: (x) => electricityColor[x.properties.__logicalLayer],
        getFilterValue: (x) => (visibility[x.properties.__logicalLayer] ? 1 : 0),
        filterRange: [1, 1],
        getLineWidth: 10,
        lineWidthUnit: 'meters',
        lineWidthMinPixels: 1,
        lineWidthMaxPixels: 10,
        lineJointRounded: true,
        lineCapRounded: true,
        updateTriggers: {
          getFilterValue: [visibility],
        },
        extensions: [new DataFilterExtension({ filterSize: 1 })],
      } as any),
  },
  {
    id: 'elec_nodes',
    type: 'MVTLayer',
    spatialType: 'vector',
    fn: ({ props, zoom }) =>
      new MVTLayer(props, {
        data: 'http://localhost:8080/data/elec_nodes.json',
        getFillColor: COLORS.electricity_high.deck,
        stroked: true,
        getLineColor: [255, 255, 255],
        lineWidthMinPixels: 1,
        ...pointRadius(zoom),
      } as any),
  },
  {
    id: 'rail_edges',
    type: 'MVTLayer',
    spatialType: 'vector',
    fn: ({ props, zoom }) =>
      new MVTLayer(props, {
        data: 'http://localhost:8080/data/rail_edges.json',
        getLineColor: COLORS.railway.deck,
        ...lineStyle(zoom),
      } as any),
  },
  {
    id: 'rail_nodes',
    type: 'MVTLayer',
    spatialType: 'vector',
    fn: ({ props, zoom }) =>
      new MVTLayer(props, {
        data: 'http://localhost:8080/data/rail_nodes.json',
        getFillColor: COLORS.railway.deck,
        stroked: true,
        getLineColor: [255, 255, 255],
        ...pointRadius(zoom),
      } as any),
  },
  {
    id: 'road_edges',
    type: 'MVTLayer',
    spatialType: 'vector',
    fn: ({ props, zoom }) =>
      new MVTLayer(props, {
        data: 'http://localhost:8080/data/road_edges.json',
        getLineColor: (x) => {
          const roadClassProp = x.properties.road_class;
          // console.log('prop', roadClassProp);
          const roadClassEnum = roadClassLookup[roadClassProp];
          // console.log('enum', roadClassEnum);
          const color = roadColor[roadClassEnum];
          // console.log(color);
          return color;
        },
        ...lineStyle(zoom),
      } as any),
  },
  {
    id: 'road_bridges',
    type: 'MVTLayer',
    spatialType: 'vector',
    fn: ({ props, zoom }) =>
      new MVTLayer(props, {
        data: 'http://localhost:8080/data/road_bridges.json',
        getFillColor: COLORS.bridges.deck,
        stroked: true,
        getLineColor: [255, 255, 255],
        ...pointRadius(zoom),
      } as any),
  },
  {
    id: 'water_potable_edges',
    type: 'MVTLayer',
    spatialType: 'vector',
    fn: ({ props, zoom }) =>
      new MVTLayer(props, {
        data: 'http://localhost:8080/data/water_potable_edges.json',
        ...lineStyle(zoom),
        getLineColor: COLORS.water_edges.deck,
      } as any),
  },
  {
    id: 'water_potable_nodes',
    type: 'MVTLayer',
    spatialType: 'vector',
    fn: ({ props, zoom }) =>
      new MVTLayer(props, {
        data: 'http://localhost:8080/data/water_potable_nodes.json',
        getFillColor: COLORS.water_abstraction.deck,
        stroked: true,
        getLineColor: [255, 255, 255],
        ...pointRadius(zoom),
      } as any),
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

  hazardDeckLayer('coastal', 1, '4x5', 2050, 'None'),
  hazardDeckLayer('coastal', 2, '4x5', 2050, 'None'),
  hazardDeckLayer('coastal', 5, '4x5', 2050, 'None'),
  hazardDeckLayer('coastal', 10, '4x5', 2050, 'None'),
  hazardDeckLayer('coastal', 50, '4x5', 2050, 'None'),
  hazardDeckLayer('coastal', 100, '4x5', 2050, 'None'),
  // {
  //   id: 'hazard',
  //   type: 'TileLayer',
  //   spatialType: 'raster',
  //   fn: ({ props, params: { floodType, returnPeriod } }) =>
  //     new TileLayer(props, {
  //       data: `http://localhost:5000/singleband/${floodType}/${returnPeriod}/raw/{z}/{x}/{y}.png?colormap=${rasterColormaps[floodType]}&stretch_range=${rasterColormapRanges[floodType]}`,
  //       refinementStrategy: 'no-overlap',
  //       renderSubLayers: (props) =>
  //         new BitmapLayer(props, {
  //           data: null,
  //           image: props.data,
  //           bounds: getBoundsForTile(props.tile),
  //           textureParameters: {
  //             [GL.TEXTURE_MIN_FILTER]: GL.NEAREST,
  //             [GL.TEXTURE_MAG_FILTER]: GL.NEAREST,
  //           },
  //         }),
  //     }),
  // },
]);

function hazardDeckLayer(hazardType, returnPeriod, rcp, epoch, confidence) {
  const id = getHazardId({ hazardType, returnPeriod, rcp, epoch, confidence }); //`hazard_${hazardType}_${returnPeriod}`;
  return {
    id,
    type: 'TileLayer',
    spatialType: 'raster',
    fn: ({ props, zoom, params: { hazardType, returnPeriod, rcp, epoch, confidence } }) =>
      new TileLayer(props, {
        data: `http://localhost:5000/singleband/${hazardType}/${returnPeriod}/${rcp}/${epoch}/${confidence}/{z}/{x}/{y}.png?colormap=${rasterColormaps[hazardType]}&stretch_range=${rasterColormapRanges[hazardType]}`,
        refinementStrategy: 'never',
        renderSubLayers: (props) =>
          new BitmapLayer(props, {
            data: null,
            image: props.data,
            bounds: getBoundsForTile(props.tile),
            textureParameters: {
              [GL.TEXTURE_MAG_FILTER]: GL.LINEAR,
              // [GL.TEXTURE_MAG_FILTER]: zoom < 12 ? GL.NEAREST : GL.NEAREST_MIPMAP_LINEAR,
            },
          }),
      }),
  };
}
