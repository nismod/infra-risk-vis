import _ from 'lodash';
import { colorCssToRgb } from '../helpers';

function makeColorConfig<K extends string>(cfg: Record<K, string>) {
  return _.mapValues(cfg, (c) => ({ css: c, deck: colorCssToRgb(c) }));
}

export const COLORS = makeColorConfig({
  electricity_high: '#eca926',
  electricity_low: '#f1d75c',
  electricity_unknown: '#aaa',

  railway: '#444',

  roads_unknown: '#b2afaa',
  roads_class_a: '#941339',
  roads_class_b: '#cb3e4e',
  roads_class_c: '#8471a8',
  roads_class_metro: '#487dbc',

  bridges: '#941339',

  water_edges: '#314386',
  water_abstraction: '#4d49bc',
});
