import React, { FC, useMemo, useState } from 'react';
import { Box, Drawer, Toolbar } from '@material-ui/core';

import { DataMap } from './map/DataMap';
import { BackgroundControl } from './controls/BackgroundControl';
import { NetworkControl } from './controls/NetworkControl';
import { useLayerSelection } from './controls/use-layer-selection';
import { ViewName, VIEWS } from './config/views';
import { BackgroundName } from './config/backgrounds';
import { LayerName } from './config/layers';
import { HazardsControl } from './controls/HazardsControl';
import { useNetworkSelection } from './controls/use-network-selection';
import { useHazardSelection } from './controls/use-hazard-selection';

interface MapViewProps {
  view: ViewName;
}

const sidebarWidth = 360;

export const MapPage: FC<MapViewProps> = ({ view }) => {
  const [background, setBackground] = useState<BackgroundName>('light');

  const viewLayerNames = useMemo<LayerName[]>(() => VIEWS[view].layers as LayerName[], [view]);
  // const layerDefinitions = useMemo(
  //   () => viewLayerNames.map((layerName) => ({ ...LAYERS[layerName], key: layerName })),
  //   [viewLayerNames],
  // );

  const { networkSelection, setNetworkSelection, networkVisibilitySet } = useNetworkSelection();
  const { hazardShow, hazardOptions, hazardParams, setSingleHazardParam, setSingleHazardShow, hazardVisibilitySet } =
    useHazardSelection();

  const visibilitySets = useMemo(
    () => [networkVisibilitySet ?? {}, hazardVisibilitySet ?? {}],
    [networkVisibilitySet, hazardVisibilitySet],
  );

  const layerSelection = useLayerSelection(viewLayerNames, visibilitySets);

  return (
    <>
      <Drawer variant="permanent">
        <Box p={3} width={sidebarWidth}>
          <Toolbar /> {/* Prevents app bar from concealing content*/}
          {view === 'overview' && (
            <>
              <NetworkControl
                networkSelection={networkSelection}
                onNetworkSelection={setNetworkSelection}
                // dataLayers={layerDefinitions}
                // layerVisibility={layerSelection}
                // onLayerVisChange={updateLayerSelection}
              />
              <HazardsControl
                hazardParams={hazardParams}
                hazardShow={hazardShow}
                hazardOptions={hazardOptions}
                onSingleHazardParam={setSingleHazardParam}
                onSingleHazardShow={setSingleHazardShow}
                // hazardSelection={hazardSelection}
                // onHazardSelection={setHazardSelection}

                // layerVisibility={layerSelection}
                // onLayerVisibilityUpdate={updateLayerSelection}
              />
            </>
          )}
          <BackgroundControl background={background} onBackgroundChange={setBackground} />
        </Box>
      </Drawer>
      <Box position="absolute" top={64} left={sidebarWidth} right={0} bottom={0}>
        <DataMap background={background} layerSelection={layerSelection} view={view} />
      </Box>
    </>
  );
};
