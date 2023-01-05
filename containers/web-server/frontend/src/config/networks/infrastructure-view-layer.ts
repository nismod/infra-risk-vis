import React from 'react';

import { InteractionTarget, VectorTarget } from '@/lib/data-map/interactions/use-interactions';
import { StyleParams, ViewLayer } from '@/lib/data-map/view-layers';
import { dataColorMap } from '@/lib/deck/props/color-map';
import { border, lineStyle, pointRadius } from '@/lib/deck/props/style';
import { fillColor, strokeColor } from '@/lib/deck/props/style';

import { assetViewLayer } from '@/config/assets/asset-view-layer';
import { assetDataAccessFunction } from '@/config/assets/data-access';
import { INFRASTRUCTURE_COLORS } from '@/config/networks/colors';
import { ExtendedAssetDetails } from '@/details/features/asset-details';
import { VectorHoverDescription } from '@/map/tooltip/VectorHoverDescription';

import { AssetViewLayerCustomFunction } from '../assets/asset-view-layer';
import { INFRASTRUCTURE_LAYER_DETAILS } from './details';
import { NETWORKS_METADATA, NetworkLayerType } from './metadata';

const roadColor = {
  road_edges_motorway: INFRASTRUCTURE_COLORS.roads_motorway.css,
  road_edges_trunk: INFRASTRUCTURE_COLORS.roads_trunk.css,
  road_edges_primary: INFRASTRUCTURE_COLORS.roads_primary.css,
  road_edges_secondary: INFRASTRUCTURE_COLORS.roads_secondary.css,
  road_edges_tertiary: INFRASTRUCTURE_COLORS.roads_unknown.css,
  road_edges_other: INFRASTRUCTURE_COLORS.roads_unknown.css,
};

function makeRoadsFn(asset_id, scaleLevel) {
  return ({ zoom, dataStyle }) => [
    strokeColor(
      dataStyle?.getColor ??
        dataColorMap(
          () => asset_id,
          (x) => roadColor[x],
        ),
    ),
    lineStyle(zoom, scaleLevel),
  ];
}
/*
function potableNodesFn({ zoom, dataStyle }) {
  return [border(), pointRadius(zoom), fillColor(dataStyle?.getColor ?? INFRASTRUCTURE_COLORS.water_supply.deck)];
}

function wastewaterNodesFn({ zoom, dataStyle }) {
  return [border(), pointRadius(zoom), fillColor(dataStyle?.getColor ?? INFRASTRUCTURE_COLORS.water_wastewater.deck)];
}
*/

/*
function electricitySourceFn({ zoom, dataStyle }) {
  return [border(), fillColor(dataStyle?.getColor ?? INFRASTRUCTURE_COLORS.electricity_high.deck), pointRadius(zoom)];
}
*/

const INFRASTRUCTURE_LAYER_FUNCTIONS: Record<NetworkLayerType, AssetViewLayerCustomFunction> = {
  power_transmission: ({ zoom, dataStyle }) => [
    strokeColor(dataStyle?.getColor ?? INFRASTRUCTURE_COLORS.electricity_high.deck),
    lineStyle(zoom, 0),
  ],
  power_distribution: ({ zoom, dataStyle }) => [
    strokeColor(dataStyle?.getColor ?? INFRASTRUCTURE_COLORS.electricity_low.deck),
    lineStyle(zoom, 1),
  ],
  // elec_edges_high: ({ zoom, dataStyle }) => [
  //   strokeColor(dataStyle?.getColor ?? INFRASTRUCTURE_COLORS.electricity_high.deck),
  //   lineStyle(zoom),
  // ],
  // elec_edges_low: ({ zoom, dataStyle }) => [
  //   strokeColor(dataStyle?.getColor ?? INFRASTRUCTURE_COLORS.electricity_low.deck),
  //   lineStyle(zoom),
  // ],
  // elec_nodes_diesel: electricitySourceFn,
  // elec_nodes_gas: electricitySourceFn,
  // elec_nodes_hydro: electricitySourceFn,
  // elec_nodes_solar: electricitySourceFn,
  // elec_nodes_wind: electricitySourceFn,
  // elec_nodes_demand: ({ zoom, dataStyle }) => [
  //   border(),
  //   fillColor(dataStyle?.getColor ?? INFRASTRUCTURE_COLORS.electricity_low.deck),
  //   pointRadius(zoom),
  // ],
  // elec_nodes_pole: ({ zoom, dataStyle }) => [
  //   border(),
  //   fillColor(dataStyle?.getColor ?? INFRASTRUCTURE_COLORS.electricity_unknown.deck),
  //   pointRadius(zoom),
  // ],
  // elec_nodes_substation: ({ zoom, dataStyle }) => [
  //   border(),
  //   fillColor(dataStyle?.getColor ?? INFRASTRUCTURE_COLORS.electricity_unknown.deck),
  //   pointRadius(zoom),
  // ],
  rail_edges: ({ zoom, dataStyle }) => [
    strokeColor(dataStyle?.getColor ?? INFRASTRUCTURE_COLORS.railway.deck),
    lineStyle(zoom),
  ],
  // rail_stations: ({ zoom, dataStyle }) => [
  //   border(),
  //   fillColor(dataStyle?.getColor ?? INFRASTRUCTURE_COLORS.railway.deck),
  //   pointRadius(zoom),
  // ],
  rail_nodes: ({ zoom, dataStyle }) => [
    border(),
    fillColor(dataStyle?.getColor ?? INFRASTRUCTURE_COLORS.railway.deck),
    pointRadius(zoom),
  ],
  road_edges_motorway: makeRoadsFn('road_edges_motorway', 0),
  road_edges_trunk: makeRoadsFn('road_edges_trunk', 0),
  road_edges_primary: makeRoadsFn('road_edges_primary', 1),
  road_edges_secondary: makeRoadsFn('road_edges_secondary', 1),
  road_edges_tertiary: makeRoadsFn('road_edges_tertiary', 2),
  /*
  road_bridges: ({ zoom, dataStyle }) => [
    border(),
    fillColor(dataStyle?.getColor ?? INFRASTRUCTURE_COLORS.bridges.deck),
    pointRadius(zoom),
  ],
  airport_runways: ({ dataStyle }) => [
    border([0, 0, 0]),
    fillColor(dataStyle?.getColor ?? INFRASTRUCTURE_COLORS.airports.deck),
  ],
  airport_terminals: ({ dataStyle }) => [
    border([0, 0, 0]),
    fillColor(dataStyle?.getColor ?? INFRASTRUCTURE_COLORS.airports.deck),
  ],
  port_areas_break: ({ dataStyle }) => [border(), fillColor(dataStyle?.getColor ?? INFRASTRUCTURE_COLORS.ports.deck)],
  port_areas_container: ({ dataStyle }) => [
    border(),
    fillColor(dataStyle?.getColor ?? INFRASTRUCTURE_COLORS.ports.deck),
  ],
  port_areas_industry: ({ dataStyle }) => [
    border(),
    fillColor(dataStyle?.getColor ?? INFRASTRUCTURE_COLORS.ports.deck),
  ],
  port_areas_silo: ({ dataStyle }) => [border(), fillColor(dataStyle?.getColor ?? INFRASTRUCTURE_COLORS.ports.deck)],
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
    strokeColor(dataStyle?.getColor ?? INFRASTRUCTURE_COLORS.water_supply.deck),
  ],
  water_irrigation_edges: ({ zoom, dataStyle }) => [
    lineStyle(zoom),
    strokeColor(dataStyle?.getColor ?? INFRASTRUCTURE_COLORS.water_irrigation.deck),
  ],
  water_irrigation_nodes: ({ zoom, dataStyle }) => [
    border(),
    pointRadius(zoom),
    fillColor(dataStyle?.getColor ?? INFRASTRUCTURE_COLORS.water_irrigation.deck),
  ],
  water_waste_sewer_gravity: ({ zoom, dataStyle }) => [
    lineStyle(zoom),
    strokeColor(dataStyle?.getColor ?? INFRASTRUCTURE_COLORS.water_wastewater.deck),
  ],
  water_waste_sewer_pressure: ({ zoom, dataStyle }) => [
    lineStyle(zoom),
    strokeColor(dataStyle?.getColor ?? INFRASTRUCTURE_COLORS.water_wastewater.deck),
  ],
  water_waste_nodes_sump: wastewaterNodesFn,
  water_waste_nodes_pump: wastewaterNodesFn,
  water_waste_nodes_relift: wastewaterNodesFn,
  water_waste_nodes_wwtp: wastewaterNodesFn,
  */
};

export function infrastructureViewLayer(
  infrastructureType: NetworkLayerType,
  styleParams: StyleParams,
): ViewLayer {
  const customFn = INFRASTRUCTURE_LAYER_FUNCTIONS[infrastructureType];
  const { label, color } = NETWORKS_METADATA[infrastructureType];

  return assetViewLayer({
    assetId: infrastructureType,
    metadata: {
      spatialType: 'vector',
      interactionGroup: 'assets',
    },
    styleParams,
    customFn,
    customDataAccessFn: assetDataAccessFunction(infrastructureType),
    renderTooltip: (hover: InteractionTarget<VectorTarget>) => {
      return React.createElement(VectorHoverDescription, {
        hoveredObject: hover,
        label,
        color,
        idValue: hover.target.feature.properties.asset_id,
      });
    },
    renderDetails: (selection: InteractionTarget<VectorTarget>) => {
      const detailsComponent = INFRASTRUCTURE_LAYER_DETAILS[infrastructureType];
      const feature = selection.target.feature;

      return React.createElement(ExtendedAssetDetails, {
        feature,
        detailsComponent,
        label,
        color,
        showRiskSection: infrastructureType !== 'rail_nodes',
      });
    },
  });
}
