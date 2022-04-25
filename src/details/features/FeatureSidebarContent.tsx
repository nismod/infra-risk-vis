import React, { FC, useEffect, useState } from 'react';
import { Box, IconButton, Typography } from '@mui/material';

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
import { ViewLayer } from 'lib/data-map/view-layers';
import { NETWORKS_METADATA } from 'config/networks/metadata';
import { ColorBox } from 'map/tooltip/content/ColorBox';
import { ApiClient } from 'lib/api-client';
import { DamagesSection } from './damages/DamagesSection';
import { Download } from '@mui/icons-material';
import { downloadFile } from 'lib/helpers';

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
  road_edges_class_a: RoadEdgeDetails,
  road_edges_class_b: RoadEdgeDetails,
  road_edges_class_c: RoadEdgeDetails,
  road_edges_metro: RoadEdgeDetails,
  road_edges_track: RoadEdgeDetails,
  road_edges_other: RoadEdgeDetails,

  water_irrigation_edges: IrrigationDetails,
  water_irrigation_nodes: IrrigationDetails,
  water_potable_nodes: WaterSupplyNodeDetails,
  water_potable_edges: WaterPipelineDetails,
  water_waste_edges: WaterPipelineDetails,
  water_waste_nodes: WastewaterNodeDetails,
};

interface FeatureSidebarContentProps {
  feature: any;
  viewLayer: ViewLayer;
}

export const FeatureSidebarContent: FC<FeatureSidebarContentProps> = ({ feature, viewLayer }) => {
  const DetailsComponent = componentMapping[viewLayer.id] ?? DefaultDetails;
  const { color, label } = NETWORKS_METADATA[viewLayer.id];

  const f = feature.properties;

  const [featureDetails] = useFeatureDetails(feature?.id);

  return (
    <Box position="relative">
      <pre id="feature_debug" style={{ display: 'none' }}>
        <code>{JSON.stringify(f, null, 2)}</code>
      </pre>
      <Typography variant="caption">
        <ColorBox color={color ?? '#333'} />
        {label}
      </Typography>
      <DetailsComponent f={f} />
      {featureDetails && (
        <>
          <IconButton
            sx={{
              position: 'absolute',
              top: 0,
              right: 0,
            }}
            title="Download CSV with feature metadata"
            onClick={() => downloadFile(makeDetailsCsv(featureDetails), 'text/csv', `feature_${feature.id}.csv`)}
          >
            <Download />{' '}
          </IconButton>
          <DamagesSection fd={featureDetails} />
        </>
      )}
    </Box>
  );
};

function makeDetailsCsv(fd) {
  return (
    'variable,value\n' +
    Object.entries(fd.properties)
      .map(([k, v]) => `${k},${v}`)
      .join('\n')
  );
}

const apiClient = new ApiClient({
  BASE: '/api',
});

function useFeatureDetails(featureId: number) {
  const [featureDetails, setFeatureDetails] = useState(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState(null);

  useEffect(() => {
    setFeatureDetails(null);
    setLoading(true);
    setError(null);
    if (featureId) {
      apiClient.features
        .featuresReadFeature({ featureId })
        .then((featureDetails) => {
          setFeatureDetails(featureDetails);
          setLoading(false);
        })
        .catch((error) => {
          setError(error);
          setLoading(false);
        });
    }
  }, [featureId]);

  return [featureDetails, loading, error];
}
