import React from 'react'
import { unique } from './helpers'

const Tooltip = (props) => {
  const source_layers = unique(props.features.map(f => f.sourceLayer));

  return (source_layers.length)? (
    <div className="tooltip-wrap">
      <div className="tooltip-body">
        {
          source_layers.map((source_layer, i) => (
            <div key={i}>
              <strong>{source_layer}</strong>
            </div>
          ))
        }
      </div>
      <span className="tooltip-triangle"></span>
    </div>
  ) : null;
}

export default Tooltip;
