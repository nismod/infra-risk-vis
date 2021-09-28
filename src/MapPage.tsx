import React, { FC, useMemo, useState } from 'react';
import { Drawer, Toolbar } from '@material-ui/core';

import { DataMap } from './map/DataMap';
import { BackgroundControl } from './controls/BackgroundControl';
import { NetworkControl } from './controls/NetworkControl';
import { useLayerSelection } from './controls/use-layer-selection';
import { ViewName, VIEWS } from './config/views';
import { BackgroundName } from './config/backgrounds';
import { LayerName, LAYERS } from './config/layers';
import { HazardsControl } from './controls/HazardsControl';

interface MapViewProps {
  view: ViewName;
}

export const MapPage: FC<MapViewProps> = ({ view }) => {
  const [background, setBackground] = useState<BackgroundName>('light');

  const viewLayerNames = useMemo<LayerName[]>(() => VIEWS[view].layers as LayerName[], [view]);
  const layerDefinitions = useMemo(
    () => viewLayerNames.map((layerName) => ({ ...LAYERS[layerName], key: layerName })),
    [viewLayerNames],
  );

  const { layerSelection, updateLayerSelection } = useLayerSelection(viewLayerNames);

  return (
    <>
      <Drawer variant="permanent">
        <Toolbar /> {/* Prevents app bar from concealing content*/}
        <div className="drawer-contents">
          {view === 'overview' && (
            <>
              <NetworkControl
                dataLayers={layerDefinitions}
                layerVisibility={layerSelection}
                onLayerVisChange={updateLayerSelection}
              />
              <HazardsControl layerVisibility={layerSelection} onLayerVisibilityUpdate={updateLayerSelection} />
            </>
          )}
          <BackgroundControl background={background} onBackgroundChange={setBackground} />
        </div>
      </Drawer>
      <div className="map-height">
        <DataMap background={background} layerSelection={layerSelection} view={view} />
      </div>
    </>
  );
};
