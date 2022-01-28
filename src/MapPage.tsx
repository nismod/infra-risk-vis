import _ from 'lodash';
import { FC, useCallback, useMemo, useState } from 'react';
import { Box, Drawer, Toolbar, Typography } from '@material-ui/core';
import ToggleButton from '@material-ui/lab/ToggleButton';
import ToggleButtonGroup from '@material-ui/lab/ToggleButtonGroup';

import { DataMap } from './map/DataMap';
import { NetworkControl } from './controls/NetworkControl';
import { ViewName } from './config/views';
import { HazardsControl } from './controls/HazardsControl';
import { useNetworkSelection } from './controls/use-network-selection';
import { useHazardSelection } from './controls/hazards/use-hazard-selection';
import { useHazardVisibility } from './controls/hazards/use-hazard-visibility';
import { useHazardParams } from './controls/hazards/use-hazard-params';
import { useDataParams } from './controls/use-data-params';
import { totalDamagesConfig } from './config/data/total-damages';

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

function getEadKey(hazard: string, rcp: string, epoch: number) {
  return `${hazard}__rcp_${rcp}__epoch_${epoch}__conf_None`;
}

function getTotalDamagesAccessor(rcp, epoch) {
  return (f) =>
    _.sum(['fluvial', 'surface', 'coastal', 'cyclone'].map((ht) => f.properties[getEadKey(ht, rcp, epoch)] ?? 0));
}

function getDamageMapStyleParams(showTotalDamages, totalDamageParams, hazardParams, selectedHazard) {
  const { rcp, epoch } = showTotalDamages ? totalDamageParams : hazardParams[selectedHazard];

  const colorField = showTotalDamages ? getTotalDamagesAccessor(rcp, epoch) : getEadKey(selectedHazard, rcp, epoch);

  return {
    colorMap: {
      colorScheme: 'damages',
      colorField,
    },
  };
}

export const MapPage: FC<MapViewProps> = ({ view }) => {
  const [mode, setMode] = useState<'input' | 'direct-damages'>('input');
  const showDamages = mode === 'direct-damages';
  const [showDamageRaster, setShowDamageRaster] = useState(true);
  const [showTotalDamages, setShowTotalDamages] = useState(false);

  const { networkSelection, setNetworkSelection, networkVisibilitySet } = useNetworkSelection();
  const forceSingleHazard = showDamages;
  const {
    hazardSelection,
    setSingleHazardShow,
    clearSelection: clearHazardSelection,
  } = useHazardSelection(forceSingleHazard);

  const { hazardOptions, hazardParams, setSingleHazardParam } = useHazardParams();

  let hazardVisibilitySet = useHazardVisibility(hazardSelection, hazardParams);
  const showHazard = !showDamages || showDamageRaster;
  hazardVisibilitySet = showHazard ? hazardVisibilitySet : {};

  const {
    params: totalDamageParams,
    options: totalDamageOptions,
    updateParam: updateTotalDamageParam,
  } = useDataParams(totalDamagesConfig);

  const handleShowTotalDamages = useCallback(
    (e, show) => {
      setShowTotalDamages(show);
      if (show) {
        clearHazardSelection();
      }
    },
    [clearHazardSelection],
  );

  const handleShowSingleHazard = useCallback(
    (hazardType, show) => {
      setSingleHazardShow(hazardType, show);
      if (show) {
        setShowTotalDamages(false);
      }
    },
    [setSingleHazardShow],
  );

  const riskMapSelectedHazard = useMemo(() => firstTrue(hazardSelection), [hazardSelection]);

  const styleParams = useMemo(() => {
    if (!showDamages || !(riskMapSelectedHazard || showTotalDamages)) return {};

    return getDamageMapStyleParams(showTotalDamages, totalDamageParams, hazardParams, riskMapSelectedHazard);
  }, [showDamages, riskMapSelectedHazard, hazardParams, showTotalDamages, totalDamageParams]);

  const layerSelection = useMemo(
    () => Object.assign({}, networkVisibilitySet ?? {}, hazardVisibilitySet ?? {}),
    [networkVisibilitySet, hazardVisibilitySet],
  );

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
                onSingleHazardShow={handleShowSingleHazard}
                showDamages={showDamages}
                showDamageRaster={showDamageRaster}
                onShowDamageRaster={setShowDamageRaster}
                showTotalDamages={showTotalDamages}
                onShowTotalDamages={handleShowTotalDamages}
                totalDamagesParams={totalDamageParams}
                totalDamagesOptions={totalDamageOptions}
                onSingleTotalDamagesParam={updateTotalDamageParam}
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
