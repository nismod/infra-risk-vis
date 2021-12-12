import React, { FC } from 'react';
import { Typography } from '@material-ui/core';

import { LayerDefinition } from '../config/layers';
import { RiskSection } from './RiskSection';
import {
  AirportDetails,
  BridgeDetails,
  DefaultDetails,
  IrrigationDetails,
  PortDetails,
  PowerDemandNodeDetails,
  PowerGenerationNodeDetails,
  PowerJunctionNodeDetails,
  PowerLineDetails,
  RailEdgeDetails,
  RailNodeDetails,
  RoadEdgeDetails,
  WastewaterNodeDetails,
  WaterPipelineDetails,
  WaterSupplyNodeDetails,
} from './detail-components';

var componentMapping = {
  airport_areas: AirportDetails,
  elec_edges_high: PowerLineDetails,
  elec_edges_low: PowerLineDetails,
  elec_nodes_junction: PowerJunctionNodeDetails,
  elec_nodes_sink: PowerDemandNodeDetails,
  elec_nodes_source: PowerGenerationNodeDetails,
  port_areas: PortDetails,
  rail_edges: RailEdgeDetails,
  rail_nodes: RailNodeDetails,
  road_bridges: BridgeDetails,
  road_edges: RoadEdgeDetails,
  water_irrigation_edges: IrrigationDetails,
  water_irrigation_nodes: IrrigationDetails,
  water_potable_nodes: WaterSupplyNodeDetails,
  water_potable_edges: WaterPipelineDetails,
  water_waste_edges: WaterPipelineDetails,
  water_waste_nodes: WastewaterNodeDetails,
};

interface FeatureSidebarContentProps {
  f: any;
  layer: LayerDefinition;
}

export const FeatureSidebarContent: FC<FeatureSidebarContentProps> = ({ f, layer }) => {
  const DetailsComponent = componentMapping[layer.id] ?? DefaultDetails;

  return (
    <>
      <pre id="feature_debug" style={{ display: 'none' }}>
        <code>{JSON.stringify(f, null, 2)}</code>
      </pre>
      <Typography variant="caption">
        <span style={{ color: layer.color ?? '#333' }}>■</span>&nbsp;{layer.label}
      </Typography>
      <DetailsComponent f={f} />
      <RiskSection f={f} />
    </>
  );
};
