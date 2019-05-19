import React from 'react'

const Tooltip = (props) => {
  const seen = {};
  return (props.features.length)? (
    <div className="flex-parent-inline flex-parent--center-cross flex-parent--column absolute bottom">
      <div className="flex-child px12 py12 bg-gray-dark color-white shadow-darken10 round txt-s w240 clip txt-truncate">
        {
          props.features.map((feature, i) => {
            const source_layer = feature.layer['source-layer'];
            if (seen[source_layer]) {
              return null;
            }
            seen[source_layer] = true;
            return (
              <div key={i}>
                <strong>{source_layer}</strong>
              </div>
            )
          })
        }
      </div>
      <span className="flex-child color-gray-dark triangle triangle--d"></span>
    </div>
  ) : null;
}

export default Tooltip;
