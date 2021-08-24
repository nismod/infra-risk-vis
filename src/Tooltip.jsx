import React from 'react';
import { titleCase } from './helpers';

const Tooltip = (props) => {
  const entries = {};

  for (const f of props.features) {
    let title = titleCase(
      f.sourceLayer.replace(/_/g, ' ').replace('edges', '').replace('nodes', '').replace('elec', 'electricity'),
    );
    let subtitle = f.properties.road_type ? '(' + f.properties.road_type + ')' : '';

    let max_value;
    let detail;

    if (!entries[f.sourceLayer] || entries[f.sourceLayer].max_value < max_value) {
      entries[f.sourceLayer] = { title, subtitle, max_value, detail };
    }
  }

  return props.features.length ? (
    <div className="tooltip-wrap">
      <div className="tooltip-body">
        {Object.values(entries).map((entry, i) => {
          return (
            <div key={i}>
              <strong>
                {entry.title} {entry.subtitle}
              </strong>
              {entry.detail}
            </div>
          );
        })}
      </div>
      <span className="tooltip-triangle"></span>
    </div>
  ) : null;
};

export default Tooltip;
