import { FC, useContext } from 'react';
import { MapContext, FlyToInterpolator } from 'react-map-gl';

import { easeCubic } from 'd3-ease';

import { useChangeEffect } from 'lib/hooks/use-change-effect';
import { ViewStateContext } from 'lib/data-map/DeckMap';

type DeckBoundingBox = [[number, number], [number, number]];

export interface BoundingBox {
  minX: number;
  minY: number;
  maxX: number;
  maxY: number;
}

interface MapBoundsFitterProps {
  boundingBox: BoundingBox;
}

export const MapBoundsFitter: FC<MapBoundsFitterProps> = ({ boundingBox }) => {
  const { viewport } = useContext(MapContext);
  const { viewState, setViewState } = useContext(ViewStateContext);

  useChangeEffect(
    () => {
      if (boundingBox != null) {
        const deckBbox: DeckBoundingBox = [
          [boundingBox.minX, boundingBox.minY],
          [boundingBox.maxX, boundingBox.maxY],
        ];
        const { latitude, longitude, zoom } = viewport.fitBounds(deckBbox, { padding: 20 });

        setViewState({
          ...viewState,
          latitude,
          longitude,
          zoom,
          transitionDuration: 1500,
          transitionInterpolator: new FlyToInterpolator() as any,
          transitionEasing: easeCubic,
        });
      }
    },
    [boundingBox, viewState, setViewState, viewport],
    [boundingBox],
  );

  return null;
};
