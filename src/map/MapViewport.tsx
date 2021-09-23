import DeckGL from 'deck.gl';
import { useMemo, useState } from 'react';
import { StaticMap } from 'react-map-gl';
import _ from 'lodash';

import { backgroundConfig, BackgroundName } from '../config/backgrounds';

const MAPBOX_KEY = 'pk.eyJ1IjoidG9tcnVzc2VsbCIsImEiOiJjaXZpMTFpdGkwMDQ1MnptcTh4ZzRzeXNsIn0.ZSvSOHSsWBQ44QNhA71M6Q';

function visible(isVisible: boolean): 'visible' | 'none' {
  return isVisible ? 'visible' : 'none';
}

function makeMapboxConfig(background: BackgroundName) {
  return {
    version: 8,
    sources: _.mapValues(backgroundConfig, (b) => b.source),
    layers: Object.values(backgroundConfig).map((b) =>
      _.merge(b.layer, { layout: { visibility: visible(background === b.id) } }),
    ),
  };
}

export const MapViewport = ({ layers, background, onHover, onClick, children }) => {
  const [viewport, setViewport] = useState({
    latitude: 18.14,
    longitude: -77.28,
    zoom: 8,
    minZoom: 3,
    maxZoom: 16,
    maxPitch: 0,
  });

  const backgroundStyle = useMemo(() => makeMapboxConfig(background), [background]);

  return (
    <DeckGL
      controller={true}
      viewState={viewport}
      onViewStateChange={({ viewState }) => setViewport(viewState)}
      layers={layers}
      pickingRadius={10}
      onHover={onHover}
      onClick={onClick}
    >
      <StaticMap mapStyle={backgroundStyle} mapboxApiAccessToken={MAPBOX_KEY} />
      {children}
    </DeckGL>
  );
};
