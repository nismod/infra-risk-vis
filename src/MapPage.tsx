import React, { FC, useMemo, useState } from 'react';
import { Box, Checkbox, Drawer, FormControlLabel, Toolbar, Typography } from '@material-ui/core';

import { DataMap } from './map/DataMap';
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

  const [showDamages, setShowDamages] = useState(false);

  const { networkSelection, setNetworkSelection, networkVisibilitySet } = useNetworkSelection();
  const { hazardShow, hazardOptions, hazardParams, setSingleHazardParam, setSingleHazardShow, hazardVisibilitySet } =
    useHazardSelection(showDamages);

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
              <FormControlLabel
                control={<Checkbox checked={showDamages} onChange={(e, checked) => setShowDamages(checked)} />}
                label="Map damages"
              ></FormControlLabel>
              <HazardsControl
                hazardParams={hazardParams}
                hazardShow={hazardShow}
                hazardOptions={hazardOptions}
                onSingleHazardParam={setSingleHazardParam}
                onSingleHazardShow={setSingleHazardShow}
                showDamages={showDamages}
                // hazardSelection={hazardSelection}
                // onHazardSelection={setHazardSelection}

                // layerVisibility={layerSelection}
                // onLayerVisibilityUpdate={updateLayerSelection}
              />
            </>
          )}
        </Box>
      </Drawer>
      <Box position="absolute" top={64} left={sidebarWidth} right={0} bottom={0}>
        <DataMap
          background={background}
          onBackground={setBackground}
          layerSelection={layerSelection}
          view={view}
          showDamages={showDamages}
        />
      </Box>
    </>
  );
};
