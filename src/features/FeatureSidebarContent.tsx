import React, { FC, useMemo } from 'react';
import { Typography } from '@material-ui/core';

import { LayerDefinition } from '../config/layers';
import { RiskSection } from './RiskSection';
import { EADChartSection } from './EADChartSection';
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
import { hazardConfig } from '../controls/use-hazard-selection';

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

function getFeatureEadData(f: any) {
  const eadData = [];
  for (const [hazardType, hazard] of Object.entries(hazardConfig)) {
    for (const rcp of hazard.paramDomains.rcp) {
      for (const epoch of hazard.paramDomains.epoch) {
        // TODO check risk data for confidence
        // for (const confidence of hazard.paramDomains.confidence) {
        const riskKey = `${hazardType}__rcp_${rcp}__epoch_${epoch}__conf_None`;
        if (f[riskKey]) {
          eadData.push({ key: riskKey, hazardType, rcp, epoch: epoch.toString(), ead: f[riskKey] });
        }
        // }
      }
    }
  }
  return eadData;
}

export const FeatureSidebarContent: FC<FeatureSidebarContentProps> = ({ f, layer }) => {
  const DetailsComponent = componentMapping[layer.id] ?? DefaultDetails;
  const eadData = useMemo(() => getFeatureEadData(f), [f]);
  return (
    <>
      <pre id="feature_debug" style={{ display: 'none' }}>
        <code>{JSON.stringify(f, null, 2)}</code>
      </pre>
      <Typography variant="caption">
        <span style={{ color: layer.color ?? '#333' }}>■</span>&nbsp;{layer.label}
      </Typography>
      <DetailsComponent f={f} />
      <RiskSection eadData={eadData} />
      <EADChartSection eadData={eadData} />
    </>
  );
};
