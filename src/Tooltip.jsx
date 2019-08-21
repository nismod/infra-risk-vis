import React from 'react'
import { unique, titleCase } from './helpers'

const Tooltip = (props) => {
  const source_layers = unique(props.features.map(f => f.sourceLayer));

  return (source_layers.length)? (
    <div className="tooltip-wrap">
      <div className="tooltip-body">
        {
          source_layers.map((source_layer, i) => (
            <div key={i}>
              <strong>{titleCase(
                source_layer.replace(/_/g," ")
                  .replace(/m(\d)/, 'm-$1')
                  .replace('4m-999m', '>4m')
                  .replace('1in', '1/')
              )}</strong>
            </div>
          ))
        }
      </div>
      <span className="tooltip-triangle"></span>
    </div>
  ) : null;
}

export default Tooltip;
