import React, { FC } from 'react';
import { VegaLite } from 'react-vega';
import { MapboxGeoJSONFeature } from 'mapbox-gl';
import { List, ListItem, ListItemText, Typography } from '@material-ui/core';

import { titleCase } from './helpers';

interface FeatureSidebarProps {
  feature: MapboxGeoJSONFeature;
}

export const FeatureSidebar: FC<FeatureSidebarProps> = ({ feature }) => {
  if (!feature) {
    return null;
  }
  const f = feature.properties;

  return (
    <div className="custom-map-control top-right selected-feature">
      <Typography variant="h6">Selected Asset</Typography>
      <pre style={{ display: 'none' }}>
        <code>{JSON.stringify(f, null, 2)}</code>
      </pre>
      <List>
        {Object.entries(f).map(([key, value]) => (
          <ListItem key={key}>
            <ListItemText
              primary={titleCase(key.replace(/_/g, ' '))}
              primaryTypographyProps={{ variant: 'caption' }}
              secondary={value === ' ' ? '-' : value || '-'}
            />
          </ListItem>
        ))}
      </List>
      {f.node_id === 'Pump_682' ? (
        <div style={{ paddingLeft: '10px' }}>
          <Typography variant="caption">Stream Flow</Typography>
          <br />
          <VegaLite
            spec={{
              width: 400,
              height: 200,
              mark: 'line',
              encoding: {
                x: { field: 'date', type: 'temporal', title: 'Date' },
                y: {
                  field: 'flow_1000m3_per_day',
                  type: 'quantitative',
                  title: 'Stream Flow (1000 m³ per day)',
                },
              },
              data: { url: 'flow_cave_river.csv' },
            }}
            // defines the actions available behind the ... menu top-right of chart
            actions={{
              export: true,
              source: false,
              compiled: false,
              editor: false,
            }}
          />
        </div>
      ) : null}
    </div>
  );
};
