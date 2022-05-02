import { ZoomIn } from '@mui/icons-material';
import { IconButton, TableCell } from '@mui/material';
import { Box } from '@mui/system';
import { ExpandableRow } from 'asset-list/ExpandableRow';
import { SortedAssetTable } from 'asset-list/SortedAssetTable';
import { ListFeature } from 'asset-list/use-sorted-features';
import { getAssetDataFormats } from 'config/assets/data-formats';
import { FeatureSidebarContent } from 'details/features/FeatureSidebarContent';
import { extendBbox } from 'lib/bounding-box';
import { colorMap } from 'lib/color-map';
import { mapFitBoundsState } from 'map/MapView';
import { ColorBox } from 'map/tooltip/content/ColorBox';
import { useCallback, useMemo } from 'react';
import { atom, useRecoilState, useRecoilValue, useSetRecoilState } from 'recoil';
import { adaptationColorSpecState, adaptationFieldSpecState, adaptationLayerSpecState } from 'state/layers/networks';

import './asset-table.css';

export const hoveredAdaptationFeatureState = atom<ListFeature>({
  key: 'hoveredAdaptationFeatureState',
  default: null,
});

export const selectedAdaptationFeatureState = atom<ListFeature>({
  key: 'selectedAdaptationFeatureState',
  default: null,
});

export const FeatureAdaptationsTable = () => {
  const layerSpec = useRecoilValue(adaptationLayerSpecState);
  const fieldSpec = useRecoilValue(adaptationFieldSpecState);
  const colorSpec = useRecoilValue(adaptationColorSpecState);

  const setHoveredFeature = useSetRecoilState(hoveredAdaptationFeatureState);
  const [selectedFeature, setSelectedFeature] = useRecoilState(selectedAdaptationFeatureState);
  const setMapFitBounds = useSetRecoilState(mapFitBoundsState);

  const handleZoomInFeature = useCallback(
    (feature: ListFeature) => feature && setMapFitBounds(extendBbox(feature.bbox, 1)),
    [setMapFitBounds],
  );

  const colorFn = useMemo(() => colorMap(colorSpec), [colorSpec]);
  const { getDataLabel, getValueFormatted } = getAssetDataFormats(fieldSpec);

  return (
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
          onMouseEnter={() => setHoveredFeature(feature)}
          onMouseLeave={() => setHoveredFeature(null)}
          expandableContent={
            <Box py={1}>
              <FeatureSidebarContent feature={feature} assetType={feature.layer} showRiskSection={false} />
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
  );
};