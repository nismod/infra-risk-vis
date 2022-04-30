import { SortedAssetTable } from 'asset-list/SortedAssetTable';
import { ListFeature } from 'asset-list/use-sorted-features';
import { extendBbox } from 'lib/bounding-box';
import { mapFitBoundsState } from 'map/MapView';
import { useCallback } from 'react';
import { atom, useRecoilState, useRecoilValue, useSetRecoilState } from 'recoil';
import { adaptationFieldSpecState, adaptationLayerSpecState } from 'state/layers/networks';

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

  const [hoveredFeature, setHoveredFeature] = useRecoilState(hoveredAdaptationFeatureState);
  const [selectedFeature, setSelectedFeature] = useRecoilState(selectedAdaptationFeatureState);
  const setMapFitBounds = useSetRecoilState(mapFitBoundsState);

  const handleSelectedFeature = useCallback(
    (feature: ListFeature) => {
      setSelectedFeature(feature);
      if (feature) {
        setMapFitBounds(extendBbox(feature.bbox, 5));
      }
    },
    [setMapFitBounds, setSelectedFeature],
  );

  return (
    <SortedAssetTable
      layerSpec={layerSpec}
      fieldSpec={fieldSpec}
      hoveredFeature={hoveredFeature}
      onHoveredFeature={setHoveredFeature}
      selectedFeature={selectedFeature}
      onSelectedFeature={handleSelectedFeature}
    />
  );
};
