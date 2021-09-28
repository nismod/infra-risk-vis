import { MVTLayer, TileLayer, BitmapLayer } from 'deck.gl';
import GL from '@luma.gl/constants';

import { COLORS } from './colors';
import { makeConfig } from '../helpers';

function makeColormap(attr: string, valueColorMap: any) {
  return (x) => valueColorMap[x.properties[attr]];
}

const lineStyle = (zoom) => ({
  getLineWidth: 10,
  lineWidthUnit: 'meters',
  lineWidthMinPixels: 1,
  lineWidthMaxPixels: 10,
  lineJointRounded: true,
  lineCapRounded: true,

  // widthScale: 2 ** (15 - zoom),
});

const pointRadius = (zoom) => ({
  getPointRadius: 20,
  pointRadiusUnit: 'meters',
  pointRadiusMinPixels: 3,
  pointRadiusMaxPixels: 20,
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
};

const rasterColormapRanges = {
  fluvial: '[0,10]',
  coastal: '[0,3.5]',
};

function getBoundsForTile(tileProps) {
  const {
    bbox: { west, south, east, north },
  } = tileProps;

  return [west, south, east, north];
}

const DECK_COLOR_TRANSPARENT = [0, 0, 0, 0];

export const DECK_LAYERS = makeConfig<any, string>([
  {
    id: 'elec_edges',
    type: 'MVTLayer',
    spatialType: 'vector',
    fn: ({ props, zoom, visibility }) =>
      new MVTLayer(props, {
        data: 'http://localhost:8080/data/elec_edges.json',
        getLineColor: (x) => {
          const lt = x.properties.line_type;

          if (lt === 'High Voltage') {
            return visibility.elec_edges_high ? COLORS.electricity_high.deck : DECK_COLOR_TRANSPARENT;
          } else if (lt === 'Low Voltage') {
            return visibility.elec_edges_low ? COLORS.electricity_low.deck : DECK_COLOR_TRANSPARENT;
          } else return COLORS.electricity_unknown;
        },
        getLineWidth: 10,
        lineWidthUnit: 'meters',
        lineWidthMinPixels: 1,
        lineWidthMaxPixels: 10,
        lineJointRounded: true,
        lineCapRounded: true,
        updateTriggers: {
          getLineColor: [visibility],
        },
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
        stroked: false,
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
        stroked: false,
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
        getLineColor: makeValueMapping((x) => x.properties.class, ['CLASS A', 'CLASS B', 'CLASS C', 'METRO'])(
          [
            COLORS.roads_class_a.deck,
            COLORS.roads_class_b.deck,
            COLORS.roads_class_c.deck,
            COLORS.roads_class_metro.deck,
          ],
          COLORS.roads_unknown.deck,
        ),
        ...lineStyle(zoom),
      } as any),
  },
  {
    id: 'bridges',
    type: 'MVTLayer',
    spatialType: 'vector',
    fn: ({ props, zoom }) =>
      new MVTLayer(props, {
        data: 'http://localhost:8080/data/bridges.json',
        getFillColor: COLORS.bridges.deck,
        stroked: false,
        ...pointRadius(zoom),
      } as any),
  },
  {
    id: 'pot_edges',
    type: 'MVTLayer',
    spatialType: 'vector',
    fn: ({ props, zoom }) =>
      new MVTLayer(props, {
        data: 'http://localhost:8080/data/pot_edges.json',
        ...lineStyle(zoom),
        getLineColor: COLORS.water_edges.deck,
      } as any),
  },
  {
    id: 'abs_nodes',
    type: 'MVTLayer',
    spatialType: 'vector',
    fn: ({ props, zoom }) =>
      new MVTLayer(props, {
        data: 'http://localhost:8080/data/abs_nodes.json',
        getFillColor: COLORS.water_abstraction.deck,
        stroked: false,
        ...pointRadius(zoom),
      } as any),
  },
  hazardDeckLayer('fluvial', 20),
  hazardDeckLayer('fluvial', 50),
  hazardDeckLayer('fluvial', 100),
  hazardDeckLayer('fluvial', 200),
  hazardDeckLayer('fluvial', 500),
  hazardDeckLayer('fluvial', 1500),
  hazardDeckLayer('coastal', 1),
  hazardDeckLayer('coastal', 2),
  hazardDeckLayer('coastal', 5),
  hazardDeckLayer('coastal', 10),
  hazardDeckLayer('coastal', 50),
  hazardDeckLayer('coastal', 100),
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

function hazardDeckLayer(hazardType, returnPeriod) {
  const id = `hazard_${hazardType}_${returnPeriod}`;
  return {
    id,
    type: 'TileLayer',
    spatialType: 'raster',
    fn: ({ props, params: { floodType, returnPeriod } }) =>
      new TileLayer(props, {
        data: `http://localhost:5000/singleband/${hazardType}/${returnPeriod}/raw/{z}/{x}/{y}.png?colormap=${rasterColormaps[hazardType]}&stretch_range=${rasterColormapRanges[floodType]}`,
        refinementStrategy: 'no-overlap',
        renderSubLayers: (props) =>
          new BitmapLayer(props, {
            data: null,
            image: props.data,
            bounds: getBoundsForTile(props.tile),
            textureParameters: {
              [GL.TEXTURE_MIN_FILTER]: GL.NEAREST,
              [GL.TEXTURE_MAG_FILTER]: GL.NEAREST,
            },
          }),
      }),
  };
}
