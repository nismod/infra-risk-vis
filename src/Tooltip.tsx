import React from 'react';
import PropTypes from 'prop-types';
import { titleCase } from './helpers';
import { MapboxGeoJSONFeature } from 'mapbox-gl';

const Tooltip = ({ features }: { features: MapboxGeoJSONFeature[] }) => {
  const entries: object = {};

  for (const f of features) {
    let title = titleCase(
      f.sourceLayer.replace(/_/g, ' ').replace('edges', '').replace('nodes', '').replace('elec', 'electricity'),
    );
    let subtitle = f.properties.road_type ? '(' + f.properties.road_type + ')' : '';

    if (!entries[f.sourceLayer]) {
      entries[f.sourceLayer] = { title, subtitle };
    }
  }

  return features.length ? (
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

Tooltip.propTypes = {
  features: PropTypes.array.isRequired,
};

export default Tooltip;
