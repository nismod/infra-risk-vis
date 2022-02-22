import { makeColorConfig } from 'lib/helpers';

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

  airports: '#ccc',
  ports: '#000',

  water_supply: '#83B4FF',
  water_irrigation: '#0091C1',
  water_wastewater: '#4d49bc',

  buildings: '#999',
});
