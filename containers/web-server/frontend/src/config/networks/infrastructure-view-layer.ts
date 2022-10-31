import { StyleParams, ViewLayer } from '@/lib/data-map/view-layers';
import { dataColorMap } from '@/lib/deck/props/color-map';
import { border, lineStyle, pointRadius } from '@/lib/deck/props/style';
import { fillColor, strokeColor } from '@/lib/deck/props/style';

import { assetViewLayer } from '@/config/assets/asset-view-layer';
import { assetDataAccessFunction } from '@/config/assets/data-access';
import { COLORS } from '@/config/colors';

import { AssetViewLayerCustomFunction } from '../assets/asset-view-layer';

const roadColor = {
  road_edges_class_a: COLORS.roads_class_a.css,
  road_edges_class_b: COLORS.roads_class_b.css,
  road_edges_class_c: COLORS.roads_class_c.css,
  road_edges_metro: COLORS.roads_metro.css,
  road_edges_track: COLORS.roads_unknown.css,
  road_edges_other: COLORS.roads_unknown.css,
};

function makeRoadsFn(asset_id) {
  return ({ zoom, dataStyle }) => [
    strokeColor(
      dataStyle?.getColor ??
        dataColorMap(
          () => asset_id,
          (x) => roadColor[x],
        ),
    ),
    lineStyle(zoom),
  ];
}

function potableNodesFn({ zoom, dataStyle }) {
  return [border(), pointRadius(zoom), fillColor(dataStyle?.getColor ?? COLORS.water_supply.deck)];
}

function wastewaterNodesFn({ zoom, dataStyle }) {
  return [border(), pointRadius(zoom), fillColor(dataStyle?.getColor ?? COLORS.water_wastewater.deck)];
}

function electricitySourceFn({ zoom, dataStyle }) {
  return [border(), fillColor(dataStyle?.getColor ?? COLORS.electricity_high.deck), pointRadius(zoom)];
}

const INFRASTRUCTURE_LAYER_FUNCTIONS: Record<string, AssetViewLayerCustomFunction> = {
  elec_edges_high: ({ zoom, dataStyle }) => [
    strokeColor(dataStyle?.getColor ?? COLORS.electricity_high.deck),
    lineStyle(zoom),
  ],
  elec_edges_low: ({ zoom, dataStyle }) => [
    strokeColor(dataStyle?.getColor ?? COLORS.electricity_low.deck),
    lineStyle(zoom),
  ],
  elec_nodes_diesel: electricitySourceFn,
  elec_nodes_gas: electricitySourceFn,
  elec_nodes_hydro: electricitySourceFn,
  elec_nodes_solar: electricitySourceFn,
  elec_nodes_wind: electricitySourceFn,
  elec_nodes_demand: ({ zoom, dataStyle }) => [
    border(),
    fillColor(dataStyle?.getColor ?? COLORS.electricity_low.deck),
    pointRadius(zoom),
  ],
  elec_nodes_pole: ({ zoom, dataStyle }) => [
    border(),
    fillColor(dataStyle?.getColor ?? COLORS.electricity_unknown.deck),
    pointRadius(zoom),
  ],
  elec_nodes_substation: ({ zoom, dataStyle }) => [
    border(),
    fillColor(dataStyle?.getColor ?? COLORS.electricity_unknown.deck),
    pointRadius(zoom),
  ],
  rail_edges: ({ zoom, dataStyle }) => [strokeColor(dataStyle?.getColor ?? COLORS.railway.deck), lineStyle(zoom)],
  rail_stations: ({ zoom, dataStyle }) => [
    border(),
    fillColor(dataStyle?.getColor ?? COLORS.railway.deck),
    pointRadius(zoom),
  ],
  rail_junctions: ({ zoom, dataStyle }) => [
    border(),
    fillColor(dataStyle?.getColor ?? COLORS.railway.deck),
    pointRadius(zoom),
  ],
  road_edges_class_a: makeRoadsFn('road_edges_class_a'),
  road_edges_class_b: makeRoadsFn('road_edges_class_b'),
  road_edges_class_c: makeRoadsFn('road_edges_class_c'),
  road_edges_metro: makeRoadsFn('road_edges_metro'),
  road_edges_track: makeRoadsFn('road_edges_track'),
  road_edges_other: makeRoadsFn('road_edges_other'),
  road_bridges: ({ zoom, dataStyle }) => [
    border(),
    fillColor(dataStyle?.getColor ?? COLORS.bridges.deck),
    pointRadius(zoom),
  ],
  airport_runways: ({ dataStyle }) => [border([0, 0, 0]), fillColor(dataStyle?.getColor ?? COLORS.airports.deck)],
  airport_terminals: ({ dataStyle }) => [border([0, 0, 0]), fillColor(dataStyle?.getColor ?? COLORS.airports.deck)],
  port_areas_break: ({ dataStyle }) => [border(), fillColor(dataStyle?.getColor ?? COLORS.ports.deck)],
  port_areas_container: ({ dataStyle }) => [border(), fillColor(dataStyle?.getColor ?? COLORS.ports.deck)],
  port_areas_industry: ({ dataStyle }) => [border(), fillColor(dataStyle?.getColor ?? COLORS.ports.deck)],
  port_areas_silo: ({ dataStyle }) => [border(), fillColor(dataStyle?.getColor ?? COLORS.ports.deck)],
  water_potable_nodes_booster: potableNodesFn,
  water_potable_nodes_catchment: potableNodesFn,
  water_potable_nodes_entombment: potableNodesFn,
  water_potable_nodes_filter: potableNodesFn,
  water_potable_nodes_intake: potableNodesFn,
  water_potable_nodes_well: potableNodesFn,
  water_potable_nodes_pump: potableNodesFn,
  water_potable_nodes_relift: potableNodesFn,
  water_potable_nodes_reservoir: potableNodesFn,
  water_potable_nodes_river_source: potableNodesFn,
  water_potable_nodes_spring: potableNodesFn,
  water_potable_nodes_tank: potableNodesFn,
  water_potable_nodes_sump: potableNodesFn,
  water_potable_nodes_tp: potableNodesFn,
  water_potable_edges: ({ zoom, dataStyle }) => [
    lineStyle(zoom),
    strokeColor(dataStyle?.getColor ?? COLORS.water_supply.deck),
  ],
  water_irrigation_edges: ({ zoom, dataStyle }) => [
    lineStyle(zoom),
    strokeColor(dataStyle?.getColor ?? COLORS.water_irrigation.deck),
  ],
  water_irrigation_nodes: ({ zoom, dataStyle }) => [
    border(),
    pointRadius(zoom),
    fillColor(dataStyle?.getColor ?? COLORS.water_irrigation.deck),
  ],
  water_waste_sewer_gravity: ({ zoom, dataStyle }) => [
    lineStyle(zoom),
    strokeColor(dataStyle?.getColor ?? COLORS.water_wastewater.deck),
  ],
  water_waste_sewer_pressure: ({ zoom, dataStyle }) => [
    lineStyle(zoom),
    strokeColor(dataStyle?.getColor ?? COLORS.water_wastewater.deck),
  ],
  water_waste_nodes_sump: wastewaterNodesFn,
  water_waste_nodes_pump: wastewaterNodesFn,
  water_waste_nodes_relift: wastewaterNodesFn,
  water_waste_nodes_wwtp: wastewaterNodesFn,
};

export function infrastructureViewLayer(assetId: string, styleParams: StyleParams): ViewLayer {
  const customFn = INFRASTRUCTURE_LAYER_FUNCTIONS[assetId];
  return assetViewLayer({
    assetId,
    metadata: {
      spatialType: 'vector',
      interactionGroup: 'assets',
    },
    styleParams,
    customFn,
    customDataAccessFn: assetDataAccessFunction(assetId),
  });
}
