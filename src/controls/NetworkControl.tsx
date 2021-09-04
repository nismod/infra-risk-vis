import React, { FC } from 'react';
import Checkbox from '@material-ui/core/Checkbox';
import FormControl from '@material-ui/core/FormControl';
import FormControlLabel from '@material-ui/core/FormControlLabel';
import FormGroup from '@material-ui/core/FormGroup';
import FormLabel from '@material-ui/core/FormLabel';
import { LayerDefinition, LayerName } from '../config/layers';

interface NetworkControlProps {
  dataLayers: (LayerDefinition & { key: LayerName })[];
  layerVisibility: Record<LayerName, boolean>;
  onLayerVisChange: (layerName: LayerName, visibility: boolean) => void;
}
const NetworkControl: FC<NetworkControlProps> = ({ dataLayers, layerVisibility, onLayerVisChange }) => (
  <FormControl component="fieldset">
    <FormLabel component="legend">Infrastructure Layers</FormLabel>
    {dataLayers.map((layerData) => {
      const layerName = layerData.key;
      const label = layerData.label;
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
                onChange={(e) => onLayerVisChange(e.target.value as LayerName, e.target.checked)}
              />
            }
            label={
              <>
                <span
                  className={layerData.linear ? 'dot line' : 'dot'}
                  style={{ backgroundColor: layerData.color }}
                ></span>
                {label}
              </>
            }
          ></FormControlLabel>
        </FormGroup>
      );
    })}
  </FormControl>
);

export default NetworkControl;