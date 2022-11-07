import { makeColorConfig } from '@/lib/helpers';

export const INFRASTRUCTURE_COLORS = makeColorConfig({
  electricity_high: '#eca926',
  electricity_low: '#f1d75c',
  electricity_unknown: '#aaa',

  railway: '#444',

  roads_unknown: '#b2afaa',
  roads_motorway: '#941339',
  roads_trunk: '#cb3e4e',
  roads_primary: '#8471a8',
  roads_secondary: '#487dbc',

  bridges: '#941339',

  airports: '#d393d3',
  ports: '#b46666',

  water_supply: '#83B4FF',
  water_irrigation: '#0091C1',
  water_wastewater: '#4d49bc',
});
