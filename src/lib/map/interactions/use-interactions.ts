import DeckGL, { Deck, PickInfo } from 'deck.gl';
import { readPixelsToArray } from '@luma.gl/core';
import _ from 'lodash';
import { useCallback, useMemo } from 'react';
import { RecoilState, useRecoilCallback, useSetRecoilState } from 'recoil';

import { ViewLayer } from 'lib/view-layers';

import { hoverState, hoverPositionState, selectionState } from './interaction-state';

export type InteractionStyle = 'vector' | 'raster';
export interface InteractionGroupConfig {
  id: string;
  type: InteractionStyle;
  pickingRadius?: number;
  pickMultiple?: boolean;
  usesAutoHighlight?: boolean;
}

export interface InteractionTarget<T> {
  interactionGroup: string;
  interactionStyle: string;

  viewLayer: string;
  logicalLayer: string;

  target: T;
}

export interface RasterTarget {
  color: [number, number, number, number];
}

function processRasterTarget(info: any): RasterTarget {
  const { bitmap, sourceLayer } = info;
  if (bitmap) {
    const pixelColor = readPixelsToArray(sourceLayer.props.image, {
      sourceX: bitmap.pixel[0],
      sourceY: bitmap.pixel[1],
      sourceWidth: 1,
      sourceHeight: 1,
      sourceType: undefined,
    });

    return pixelColor[3]
      ? {
          color: pixelColor,
        }
      : null;
  }
}
export interface VectorTarget {
  feature: any;
}

function processVectorTarget(info: PickInfo<any>): VectorTarget {
  const { object } = info;

  return object
    ? {
        feature: object,
      }
    : null;
}

function processTargetByType(type: InteractionStyle, info: PickInfo<any>) {
  return type === 'raster' ? processRasterTarget(info) : processVectorTarget(info);
}

function processPickedObject(
  info: PickInfo<any>,
  type: InteractionStyle,
  groupName: string,
  viewLayerLookup: Record<string, ViewLayer>,
) {
  const layerId = info.layer.id;
  const target = processTargetByType(type, info);

  return (
    target && {
      interactionGroup: groupName,
      interactionStyle: type,
      viewLayer: layerId,
      logicalLayer: viewLayerLookup[layerId].getLogicalLayer?.({ deckLayerId: layerId, target }) ?? layerId,
      target,
    }
  );
}

function useSetInteractionGroupState(
  state: (group: string) => RecoilState<InteractionTarget<any> | InteractionTarget<any>[]>,
) {
  return useRecoilCallback(({ set }) => {
    return (groupName: string, value: InteractionTarget<any> | InteractionTarget<any>[]) => {
      set(state(groupName), value);
    };
  });
}

export function useInteractions(viewLayers: ViewLayer[], interactionGroups: InteractionGroupConfig[]) {
  const setHoverXY = useSetRecoilState(hoverPositionState);

  const setInteractionGroupHover = useSetInteractionGroupState(hoverState);
  const setInteractionGroupSelection = useSetInteractionGroupState(selectionState);

  const interactionGroupLookup = useMemo(() => _.keyBy(interactionGroups, 'id'), [interactionGroups]);

  const primaryGroup = interactionGroups[0].id;
  const primaryGroupPickingRadius = interactionGroupLookup[primaryGroup].pickingRadius;

  const interactiveLayers = useMemo(() => viewLayers.filter((x) => x.interactionGroup), [viewLayers]);
  const viewLayerLookup = useMemo(() => _.keyBy(interactiveLayers, (layer) => layer.id), [interactiveLayers]);
  const activeGroups = useMemo(
    () => _.groupBy(interactiveLayers, (layer) => layer.interactionGroup),
    [interactiveLayers],
  );

  const onHover = useCallback(
    (info: any, deck: Deck) => {
      const { x, y } = info;

      for (const [groupName, layers] of Object.entries(activeGroups)) {
        const layerIds = layers.map((layer) => layer.id);
        const interactionGroup = interactionGroupLookup[groupName];
        const { type, pickingRadius: radius, pickMultiple } = interactionGroup;

        const pickingParams = { x, y, layerIds, radius };

        if (pickMultiple) {
          const pickedObjects: PickInfo<any>[] = deck.pickMultipleObjects(pickingParams);
          const interactionTargets: InteractionTarget<any>[] = pickedObjects
            .map((info) => processPickedObject(info, type, groupName, viewLayerLookup))
            .filter(Boolean);

          setInteractionGroupHover(groupName, interactionTargets);
        } else {
          const info: PickInfo<any> = deck.pickObject(pickingParams);
          let interactionTarget: InteractionTarget<any> =
            info && processPickedObject(info, type, groupName, viewLayerLookup);

          setInteractionGroupHover(groupName, interactionTarget);
        }
      }

      setHoverXY([x, y]);
    },
    [activeGroups, interactionGroupLookup, setHoverXY, setInteractionGroupHover, viewLayerLookup],
  );

  const onClick = useCallback(
    (info: any, deck: DeckGL) => {
      const { x, y } = info;
      for (const [groupName, layers] of Object.entries(activeGroups)) {
        const layerIds = layers.map((layer) => layer.id);

        const interactionGroup = interactionGroupLookup[groupName];
        const { type, pickingRadius: radius } = interactionGroup;

        // currently only supports selecting vector features
        if (interactionGroup.type === 'vector') {
          const info = deck.pickObject({ x, y, layerIds, radius });
          let selectionTarget = info && processPickedObject(info, type, groupName, viewLayerLookup);

          setInteractionGroupSelection(groupName, selectionTarget);
        }
      }
    },
    [activeGroups, interactionGroupLookup, setInteractionGroupSelection, viewLayerLookup],
  );

  /**
   * Interaction groups which should be rendered during the hover picking pass
   */
  const hoverPassGroups = useMemo(
    () =>
      new Set(
        _.filter(interactionGroups, (group) => group.id === primaryGroup || group.usesAutoHighlight).map(
          (group) => group.id,
        ),
      ),
    [interactionGroups, primaryGroup],
  );

  const layerFilter = ({ layer, renderPass }) => {
    if (renderPass === 'picking:hover') {
      const viewLayer = viewLayerLookup[layer.id];
      return hoverPassGroups.has(viewLayer.interactionGroup);
    }
    return true;
  };

  return {
    onHover,
    onClick,
    layerFilter,
    pickingRadius: primaryGroupPickingRadius,
  };
}
