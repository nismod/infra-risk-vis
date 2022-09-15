import { makeColorConfig } from 'lib/helpers';

export const COLORS = makeColorConfig({
  electricity_high: '#eca926',
  electricity_low: '#f1d75c',
  electricity_unknown: '#aaa',

  railway: '#444',
  rail_past: '#b2afaa',
  rail_future: '#e2dfda',

  roads_unknown: '#b2afaa',
  roads_motorway: '#741329',
  roads_trunk: '#941339',
  roads_primary: '#cb3e4e',
  roads_secondary: '#8471a8',
  roads_tertiary: '#487dbc',

  bridges: '#941339',

  airports: '#d393d3',
  ports: '#b46666',

  water_supply: '#83B4FF',
  water_irrigation: '#0091C1',
  water_wastewater: '#4d49bc',

  buildings_commercial: '#f2808c',
  buildings_industrial: '#cb97f4',
  buildings_institutional: '#808cf2',
  buildings_mixed: '#f09e69',
  buildings_other: '#bfb4c2',
  buildings_recreation: '#95e78b',
  buildings_residential: '#f2e680',
  buildings_resort: '#f9b2ea',
  buildings_unknown: '#dfe4de',

  regions_no_data: '#ccc',
});
