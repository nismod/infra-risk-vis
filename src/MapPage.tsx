import React, { FC, useMemo, useState } from 'react';
import { Box, Drawer, Toolbar, Typography } from '@material-ui/core';
import ToggleButton from '@material-ui/lab/ToggleButton';
import ToggleButtonGroup from '@material-ui/lab/ToggleButtonGroup';

import { DataMap } from './map/DataMap';
import { NetworkControl } from './controls/NetworkControl';
import { useLayerSelection } from './controls/use-layer-selection';
import { ViewName, VIEWS } from './config/views';
import { LayerName } from './config/layers';
import { HazardsControl } from './controls/HazardsControl';
import { useNetworkSelection } from './controls/use-network-selection';
import { useHazardSelection } from './controls/hazards/use-hazard-selection';
import { useHazardVisibility } from './controls/hazards/use-hazard-visibility';
import { useHazardParams } from './controls/hazards/use-hazard-params';

interface MapViewProps {
  view: ViewName;
}

const sidebarWidth = 360;

function firstTrue(object) {
  const trueKeys = Object.entries(object)
    .filter(([key, value]) => value)
    .map(([key]) => key);
  return trueKeys[0];
}

export const MapPage: FC<MapViewProps> = ({ view }) => {
  const viewLayerNames = useMemo<LayerName[]>(() => VIEWS[view].layers as LayerName[], [view]);
  // const layerDefinitions = useMemo(
  //   () => viewLayerNames.map((layerName) => ({ ...LAYERS[layerName], key: layerName })),
  //   [viewLayerNames],
  // );

  const [mode, setMode] = useState<'input' | 'direct-damages'>('input');
  const showDamages = mode === 'direct-damages';
  const [showDamageRaster, setShowDamageRaster] = useState(true);

  const { networkSelection, setNetworkSelection, networkVisibilitySet } = useNetworkSelection();
  const forceSingleHazard = showDamages;
  const { hazardSelection, setSingleHazardShow } = useHazardSelection(forceSingleHazard);

  const { hazardOptions, hazardParams, setSingleHazardParam } = useHazardParams();

  const hazardVisibilitySet = useHazardVisibility(
    hazardSelection,
    hazardParams,
    forceSingleHazard && !showDamageRaster,
  );

  const riskMapSelectedHazard = useMemo(
    () => (showDamages ? firstTrue(hazardSelection) : null),
    [showDamages, hazardSelection],
  );

  const styleParams = useMemo(() => {
    if (!riskMapSelectedHazard) return {};

    const { rcp, epoch } = hazardParams[riskMapSelectedHazard];

    return {
      colorMap: {
        colorScheme: 'damages',
        colorField: `${riskMapSelectedHazard}__rcp_${rcp}__epoch_${epoch}__conf_None`,
      },
    };
  }, [riskMapSelectedHazard, hazardParams]);

  const visibilitySets = useMemo(
    () => [networkVisibilitySet ?? {}, hazardVisibilitySet ?? {}],
    [networkVisibilitySet, hazardVisibilitySet],
  );

  const layerSelection = useLayerSelection(viewLayerNames, visibilitySets);

  return (
    <>
      <Drawer variant="permanent">
        <Box pt={2} p={4} width={sidebarWidth} boxSizing="border-box">
          <Toolbar /> {/* Prevents app bar from concealing content*/}
          {view === 'overview' && (
            <>
              <NetworkControl networkSelection={networkSelection} onNetworkSelection={setNetworkSelection} />
              <Box mb={1}>
                <Typography variant="h6">View Mode</Typography>
              </Box>
              <ToggleButtonGroup
                size="small"
                exclusive
                value={mode}
                onChange={(e, value) => setMode(value)}
                style={{ display: 'flex', flexDirection: 'row', justifyItems: 'stretch' }}
              >
                <ToggleButton style={{ width: '50%' }} value="input">
                  Input Data
                </ToggleButton>
                <ToggleButton style={{ width: '50%' }} value="direct-damages">
                  Direct Damages
                </ToggleButton>
              </ToggleButtonGroup>
              <HazardsControl
                hazardParams={hazardParams}
                hazardShow={hazardSelection}
                hazardOptions={hazardOptions}
                onSingleHazardParam={setSingleHazardParam}
                onSingleHazardShow={setSingleHazardShow}
                showDamages={showDamages}
                showDamageRaster={showDamageRaster}
                onShowDamageRaster={setShowDamageRaster}
                forceSingle={forceSingleHazard}
              />
            </>
          )}
        </Box>
      </Drawer>
      <Box position="absolute" top={64} left={sidebarWidth} right={0} bottom={0}>
        <DataMap layerSelection={layerSelection} styleParams={styleParams} view={view} />
      </Box>
    </>
  );
};
