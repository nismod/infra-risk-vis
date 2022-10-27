import { ViewLayer } from '@/lib/data-map/view-layers';
import { selectableMvtLayer } from '@/lib/deck/layers/selectable-mvt-layer';
import { border, fillColor, pointRadius, setAlpha } from '@/lib/deck/props/style';

import { SOURCES } from '../sources';
import { PROTECTED_AREA_COLORS, ProtectedAreaType } from './metadata';

export function protectedAreaViewLayer(type: ProtectedAreaType): ViewLayer {
  const color = PROTECTED_AREA_COLORS[type];

  const id = `wdpa_${type}`;
  return {
    id,
    spatialType: 'vector',
    interactionGroup: 'wdpa',
    params: {
      type,
    },
    fn: ({ deckProps, zoom, selection }) => {
      const selectionId = selection?.target.feature.properties['WDPA_PID'];

      return [
        selectableMvtLayer(
          {
            selectionOptions: {
              selectedFeatureId: selection?.deckLayerId === `${id}@polygons` ? selectionId : undefined,
              uniqueIdProperty: 'WDPA_PID',
              selectionFillColor: [0, 0, 0, 0],
            },
          },
          deckProps,
          {
            id: `${id}@polygons`,
            data: SOURCES.vector.getUrl(`${deckProps.id}_polygons`),
            uniqueIdProperty: 'WDPA_PID',
            refinementStrategy: 'no-overlap',
            filled: true,
          },
          fillColor(setAlpha(color.deck, 100)),
          {
            highlightColor: [255, 255, 255, 100],
          },
        ),
        selectableMvtLayer(
          {
            selectionOptions: {
              selectedFeatureId: selection?.deckLayerId === `${id}@points` ? selectionId : undefined,
              uniqueIdProperty: 'WDPA_PID',
            },
          },
          deckProps,
          {
            id: `${id}@points`,
            data: SOURCES.vector.getUrl(`${deckProps.id}_points`),
            uniqueIdProperty: 'WDPA_PID',
            filled: true,
          },
          pointRadius(zoom),
          fillColor(color.deck),
          border(),
        ),
      ];
    },
  };
}