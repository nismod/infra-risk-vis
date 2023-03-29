import { AttributionControl, NavigationControl, ScaleControl } from 'react-map-gl';

import { withProps } from 'lib/react/with-props';

export const MapHudNavigationControl = withProps(
  NavigationControl,
  {
    style: {
      position: 'static',
    },
  },
  'MapHudNavigationControl',
);

export const MapHudAttributionControl = withProps(
  AttributionControl,
  {
    style: {
      position: 'static',
    },
  },
  'MapHudAttributionControl',
);

export const MapHudScaleControl = withProps(
  ScaleControl,
  {
    style: {
      position: 'static',
    },
  },
  'MapHudScaleControl',
);
