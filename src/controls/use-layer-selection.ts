import { useCallback, useState } from 'react';
import { LayerName } from '../config/layers';

export function useLayerSelection(layers: LayerName[]) {
  const [layerSelection, setLayerSelection] = useState<Record<LayerName, boolean>>(
    Object.fromEntries(layers.map((l) => [l, true])) as Record<LayerName, boolean>,
  );

  const updateLayerSelection = useCallback(
    (layerName: LayerName, selected: boolean) => {
      setLayerSelection({ ...layerSelection, [layerName]: selected });
    },
    [layerSelection],
  );

  const selectSingleLayer = useCallback(
    (layerName: LayerName) => {
      setLayerSelection({
        ...(Object.fromEntries(layers.map((l) => [l, false])) as Record<LayerName, boolean>),
        [layerName]: true,
      });
    },
    [layers],
  );

  return {
    layerSelection,
    updateLayerSelection,
    selectSingleLayer,
  };
}
