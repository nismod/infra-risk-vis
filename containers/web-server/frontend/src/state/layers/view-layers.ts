import bboxPolygon from '@turf/bbox-polygon';
import { selector } from 'recoil';

import { extendBbox } from '@/lib/bounding-box';
import { ViewLayer, viewOnlyLayer } from '@/lib/data-map/view-layers';
import { boundingBoxLayer } from '@/lib/deck/layers/bounding-box-layer';
import { truthyKeys } from '@/lib/helpers';
import { ConfigTree } from '@/lib/nested-config/config-tree';

import { buildingsViewLayer } from '@/config/buildings/buildings-view-layer';
import { labelsLayer } from '@/config/deck-layers/labels-layer';
import { hoveredAdaptationFeatureState } from '@/details/adaptations/FeatureAdaptationsTable';
import { showLabelsState } from '@/map/layers/layers-state';
import { buildingSelectionState } from '@/state/buildings';
import { isRetinaState } from '@/state/is-retina';
import { sectionVisibilityState } from '@/state/sections';

import { hazardLayerState } from './hazards';
import { networkLayersState } from './networks';
import { populationLayerState } from './population';

const buildingLayersState = selector<ViewLayer[]>({
  key: 'buildingLayersState',
  get: ({ get }) =>
    get(sectionVisibilityState('buildings'))
      ? truthyKeys(get(buildingSelectionState)).map((buildingType) => buildingsViewLayer(buildingType))
      : [],
});

export const featureBoundingBoxLayerState = selector<ViewLayer>({
  key: 'featureBoundingBoxLayerState',
  get: ({ get }) => {
    const hoveredAdaptationFeature = get(hoveredAdaptationFeatureState);

    if (!hoveredAdaptationFeature) return null;

    const geom = bboxPolygon(extendBbox(hoveredAdaptationFeature.bbox, 5));

    return viewOnlyLayer(`feature-bounding-box-${hoveredAdaptationFeature.id}`, ({ deckProps }) =>
      boundingBoxLayer({ bboxGeom: geom }, deckProps),
    );
  },
});

export const viewLayersState = selector<ConfigTree<ViewLayer>>({
  key: 'viewLayersState',
  get: ({ get }) => {
    const showLabels = get(showLabelsState);
    const isRetina = get(isRetinaState);

    return [
      get(populationLayerState),

      // hazard data layers
      get(hazardLayerState),

      get(buildingLayersState),

      // network data layers
      get(networkLayersState),

      get(featureBoundingBoxLayerState),

      showLabels && [
        // basemap labels
        viewOnlyLayer('labels', () => labelsLayer(isRetina)),
      ],

      /**
       * CAUTION: for some reason, vector layers put here are obscured by the 'labels' semi-transparent raster layer
       */
    ];
  },
});
