import * as d3 from 'd3-scale';
import { colorCssToRgb } from '../../helpers';
import { VECTOR_COLOR_MAPS } from '../color-maps';

export const mergeUpdateTriggers = (...propsArray) => {
  const updateTriggers = {};
  for (const props of propsArray) {
    const triggers = props.updateTriggers;
    if (triggers) {
      for (const [key, value] of Object.entries(triggers)) {
        updateTriggers[key] = value;
      }
    }
  }
  return { updateTriggers };
};

export const lineStyle = (zoom) => ({
  getLineWidth: 15,
  lineWidthUnit: 'meters',
  lineWidthMinPixels: 1,
  lineWidthMaxPixels: 5,
  lineJointRounded: true,
  lineCapRounded: true,

  // widthScale: 2 ** (15 - zoom),
});

export const pointRadius = (zoom) => ({
  getPointRadius: 20,
  pointRadiusUnit: 'meters',
  pointRadiusMinPixels: 3,
  pointRadiusMaxPixels: 10,
  // radiusScale: 2 ** (15 - zoom),
});

export interface ColorMapDefinition {
  colorScheme: string;
  colorField: string;
}

function makeColorMap(definition: ColorMapDefinition) {
  const { colorScheme, colorField } = definition;
  const { scale, range, empty } = VECTOR_COLOR_MAPS[colorScheme];

  const scaleFn = d3.scaleSequential(range, scale);

  return (f) => {
    const value = f.properties[colorField];
    return colorCssToRgb(value == null || value === 0 ? empty : scaleFn(value));
  };
}

export function vectorColor(type: 'fill' | 'stroke', defaultValue, styleParams) {
  const prop = styleParams?.colorMap ? makeColorMap(styleParams.colorMap) : defaultValue;

  if (type === 'fill') return { getFillColor: prop, updateTriggers: { getFillColor: [styleParams] } };
  else if (type === 'stroke') return { getLineColor: prop, updateTriggers: { getLineColor: [styleParams] } };
}

export function border(color = [255, 255, 255]) {
  return {
    stroked: true,
    getLineColor: color,
    lineWidthMinPixels: 1,
  };
}
