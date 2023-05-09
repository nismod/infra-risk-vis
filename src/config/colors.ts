import { makeColorConfig } from 'lib/helpers';

export const COLORS = makeColorConfig({
  electricity_high: '#eca926',
  electricity_low: '#f1d75c',

  elec_nodes_diesel: '#fdb1ab',
  elec_nodes_gas: '#e0c8e1',
  elec_nodes_hydro: '#b5cae0',
  elec_nodes_solar: '#ffd6a3',
  elec_nodes_wind: '#cee8c2',

  electricity_demand: '#ff8c00',
  electricity_unknown: '#aaaaaa',

  railway: '#444',

  roads_unknown: '#b2afaa',
  roads_class_a: '#941339',
  roads_class_b: '#cb3e4e',
  roads_class_c: '#8471a8',
  roads_metro: '#487dbc',

  bridges: '#941339',

  airport_runways: '#d393d3',
  airport_terminals: '#b393d3',

  port_areas_break: '#b46666',
  port_areas_container: '#b4667a',
  port_areas_industry: '#b47a66',
  port_areas_silo: '#b48e66',

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
