import { ZoomIn, ZoomOut } from '@mui/icons-material';
import { IconButton, TableCell } from '@mui/material';
import { Box } from '@mui/system';
import { useCallback, useMemo } from 'react';
import { atom, useRecoilState, useRecoilValue, useSetRecoilState } from 'recoil';

import { BoundingBox, extendBbox } from '@/lib/bounding-box';
import { colorMap } from '@/lib/color-map';
import { ColorBox } from '@/lib/ui/data-display/ColorBox';

import { ExpandableRow } from '@/asset-list/ExpandableRow';
import { SortedAssetTable } from '@/asset-list/SortedAssetTable';
import { ListFeature } from '@/asset-list/use-sorted-features';
import { getAssetDataFormats } from '@/config/assets/data-formats';
import { mapFitBoundsState } from '@/map/MapView';
import {
  adaptationColorSpecState,
  adaptationFieldSpecState,
  adaptationLayerSpecState,
} from '@/state/data-selection/networks/adaptations';

import './asset-table.css';

import { boundedFeatureState } from '@/state/layers/ui-layers/feature-bbox';

export const selectedAdaptationFeatureState = atom<ListFeature>({
  key: 'selectedAdaptationFeatureState',
  default: null,
});

const GLOBAL_BBOX: BoundingBox = [-175, -85, 179, 85];

export const FeatureAdaptationsTable = () => {
  const layerSpec = useRecoilValue(adaptationLayerSpecState);
  const fieldSpec = useRecoilValue(adaptationFieldSpecState);
  const colorSpec = useRecoilValue(adaptationColorSpecState);

  const setHoveredFeature = useSetRecoilState(boundedFeatureState);
  const [selectedFeature, setSelectedFeature] = useRecoilState(selectedAdaptationFeatureState);
  const setMapFitBounds = useSetRecoilState(mapFitBoundsState);

  const handleZoomInFeature = useCallback(
    (feature: ListFeature) => feature && setMapFitBounds(extendBbox(feature.bbox, 1)),
    [setMapFitBounds],
  );

  const handleZoomOutGlobal = useCallback(
    () => setMapFitBounds([...GLOBAL_BBOX]),
    [setMapFitBounds],
  );

  const colorFn = useMemo(() => colorMap(colorSpec), [colorSpec]);
  const { getDataLabel, getValueFormatted } = getAssetDataFormats(fieldSpec);

  return (
    <>
      <Box position="absolute" top={0} right={25} zIndex={1000}>
        <IconButton onClick={handleZoomOutGlobal} title="Zoom out to globe">
          <ZoomOut />
        </IconButton>
      </Box>
      <SortedAssetTable
        layerSpec={layerSpec}
        fieldSpec={fieldSpec}
        header={
          <>
            <TableCell width={10}>#</TableCell>
            <TableCell>{getDataLabel(fieldSpec)}</TableCell>
            <TableCell width={10}> </TableCell>
          </>
        }
        renderRow={(feature, localIndex, globalIndex) => (
          <ExpandableRow
            key={feature.string_id}
            expanded={feature === selectedFeature}
            onExpandedChange={(expanded) => setSelectedFeature(expanded ? feature : null)}
            onMouseEnter={() => setHoveredFeature({ id: `${feature.id}`, bbox: feature.bbox })}
            onMouseLeave={() => setHoveredFeature(null)}
            expandableContent={
              <Box py={1}>
                {/* TODO: render details based on view layer */}
                {/* <FeatureSidebarContent feature={feature} assetType={feature.layer} showRiskSection={false} /> */}
              </Box>
            }
          >
            <TableCell>{globalIndex + 1}</TableCell>
            <TableCell>
              <ColorBox color={colorFn(feature.value)} />
              {getValueFormatted(feature.value, fieldSpec)}
            </TableCell>
            <TableCell>
              <IconButton
                title="Zoom in to asset"
                className="row-hovered-visible"
                size="small"
                sx={{
                  padding: 0,
                }}
                onClick={(e) => {
                  handleZoomInFeature(feature);
                  e.stopPropagation();
                }}
              >
                <ZoomIn />
              </IconButton>
            </TableCell>
          </ExpandableRow>
        )}
      />
    </>
  );
};
