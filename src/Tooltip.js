import React from 'react'

const Tooltip = (props) => {
  return (props.features.length)? (
    <div className="flex-parent-inline flex-parent--center-cross flex-parent--column absolute bottom">
      <div className="flex-child px12 py12 bg-gray-dark color-white shadow-darken10 round txt-s w240 clip txt-truncate">
        {
          props.features.map((feature, i) => (
            <div key={i}>
              <strong>{feature.layer['source-layer']}</strong>
            </div>
          ))
        }
      </div>
      <span className="flex-child color-gray-dark triangle triangle--d"></span>
    </div>
  ) : null;
}

export default Tooltip;
