import React, { FC } from 'react';
import { MapboxGeoJSONFeature } from 'mapbox-gl';
import { Marker } from 'react-map-gl';

import { titleCase } from '../helpers';

interface MapTooltipProps {
  features: MapboxGeoJSONFeature[];
  tooltipLngLat: [number, number];
}

export const MapTooltip: FC<MapTooltipProps> = ({ features, tooltipLngLat }) => {
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

  return features.length && tooltipLngLat ? (
    <Marker longitude={tooltipLngLat[0]} latitude={tooltipLngLat[1]} offsetLeft={-150}>
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
    </Marker>
  ) : null;
};
