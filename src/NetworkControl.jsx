import React, { Fragment } from 'react';
import PropTypes from 'prop-types';

const NetworkControl = (props) => (
  <Fragment>
    <h4 className="h5">Networks</h4>
    {
      props.dataLayers.map(layer_data => {
        const layer = layer_data.key;
        const label = layer_data.label
        return (
          <div className="form-check" key={'toggleLayer' + layer} >
            <input className="form-check-input"
              type="checkbox"
              data-layer={layer}
              defaultChecked={true}
              id={'toggleLayerCheckbox' + layer}
              onClick={props.onLayerVisChange}/>
            <span
              className={layer_data.linear? 'dot line': 'dot'}
              style={{backgroundColor: layer_data.color}}></span>
            <label className="form-check-label" htmlFor={'toggleLayerCheckbox' + layer}>
              {label}
            </label>
          </div>
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
