import React, { FC } from 'react';
import { Checkbox, FormControl, FormControlLabel, FormGroup, FormLabel } from '@material-ui/core';

import { LayerDefinition, LayerName, LAYERS } from '../config/layers';

interface NetworkControlProps {
  dataLayers: (LayerDefinition & { key: LayerName })[];
  layerVisibility: Record<LayerName, boolean>;
  onLayerVisChange: (visUpdate: Record<string, boolean>) => void; //TODO change record key type to LayerName
}

const networkLayers = [
  'elec_edges_high',
  'elec_edges_low',
  'elec_nodes',
  'rail_edges',
  'rail_nodes',
  'road_edges',
  'bridges',
  'pot_edges',
  'abs_nodes',
];

export const NetworkControl: FC<NetworkControlProps> = ({ dataLayers, layerVisibility, onLayerVisChange }) => (
  <FormControl component="fieldset">
    <FormLabel component="legend">Infrastructure Layers</FormLabel>
    {networkLayers.map((layerName) => {
      const { label, type, color } = LAYERS[layerName];
      const checked = layerVisibility[layerName];
      return (
        <FormGroup row key={'toggleLayer' + layerName}>
          <FormControlLabel
            control={
              <Checkbox
                data-layer={layerName}
                color="primary"
                checked={checked}
                value={layerName}
                name={'toggleLayerCheckbox' + layerName}
                onChange={(e) => onLayerVisChange({ [e.target.value as LayerName]: e.target.checked })}
              />
            }
            label={
              <>
                <span className={type === 'line' ? 'dot line' : 'dot'} style={{ backgroundColor: color }}></span>
                {label}
              </>
            }
          ></FormControlLabel>
        </FormGroup>
      );
    })}
  </FormControl>
);
