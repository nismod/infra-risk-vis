import bboxPolygon from '@turf/bbox-polygon';
import { atom, selector } from 'recoil';

import { BoundingBox, extendBbox } from '@/lib/bounding-box';
import { ViewLayer, viewOnlyLayer } from '@/lib/data-map/view-layers';
import { boundingBoxLayer } from '@/lib/deck/layers/bounding-box-layer';

export interface FeatureWithBbox {
  bbox: BoundingBox;
  id: string;
}

export const boundedFeatureState = atom<FeatureWithBbox>({
  key: 'boundedFeatureBboxState',
  default: null,
});

export const featureBoundingBoxLayerState = selector<ViewLayer>({
  key: 'featureBoundingBoxLayerState',
  get: ({ get }) => {
    const { id, bbox } = get(boundedFeatureState) ?? {};

    if (!bbox) return null;

    const geom = bboxPolygon(extendBbox(bbox, 5));

    return viewOnlyLayer(`feature-bounding-box-${id}`, ({ deckProps }) =>
      boundingBoxLayer({ bboxGeom: geom }, deckProps),
    );
  },
});
