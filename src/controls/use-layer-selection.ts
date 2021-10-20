import { useMemo } from 'react';
// import { useCallback, useMemo, useState } from 'react';
import { LayerName } from '../config/layers';

export type VisibilitySet = Record<string, boolean>;

export function useLayerSelection(layers: LayerName[], visibilitySets: VisibilitySet[]) {
  //Object.fromEntries(layers.map((l) => [l, false]));

  const layerSelection = useMemo(
    () => {
      const baseLayerSelection = {};
      return Object.assign(baseLayerSelection, ...visibilitySets);
    },
    [visibilitySets],
  );

  // const updateLayerSelection = useCallback(
  //   (selectionUpdate: Record<string, boolean>) => {
  //     setLayerSelection({ ...layerSelection, ...selectionUpdate });
  //   },
  //   [layerSelection],
  // );

  // const selectSingleLayer = useCallback(
  //   (layerName: LayerName) => {
  //     setLayerSelection({
  //       ...(Object.fromEntries(layers.map((l) => [l, false])) as Record<LayerName, boolean>),
  //       [layerName]: true,
  //     });
  //   },
  //   [layers],
  // );

  // return {
  //   layerSelection,
  //   updateLayerSelection,
  //   selectSingleLayer,
  // };

  return layerSelection;
}
