import React, { Fragment } from 'react';
import PropTypes from 'prop-types';
import Checkbox from '@material-ui/core/Checkbox';
import FormControlLabel from '@material-ui/core/FormControlLabel';
import FormGroup from '@material-ui/core/FormGroup';

const NetworkControl = (props) => (
  <Fragment>
    {
      props.dataLayers.map(layer_data => {
        const layer = layer_data.key;
        const label = layer_data.label;
        const checked = props.layerVisibility[layer];
        return (
          <FormGroup row key={'toggleLayer' + layer}>
          <FormControlLabel
            control={
              <Checkbox
                data-layer={layer}
                color="primary"
                checked={checked}
                value={layer}
                name={'toggleLayerCheckbox' + layer}
                onChange={props.onLayerVisChange} />
            }
            label={
              <Fragment>
                <span
                  className={layer_data.linear? 'dot line': 'dot'}
                  style={{backgroundColor: layer_data.color}}></span>
                {label}
              </Fragment>
            }
            >
          </FormControlLabel>
          </FormGroup>
        )
      })
    }
  </Fragment>
)

NetworkControl.propTypes = {
  onLayerVisChange: PropTypes.func,
  dataLayers: PropTypes.arrayOf(PropTypes.shape({
    key: PropTypes.string,
    label: PropTypes.string
  }))
}

export default NetworkControl;
