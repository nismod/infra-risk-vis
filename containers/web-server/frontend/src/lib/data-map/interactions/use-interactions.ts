import { readPixelsToArray } from '@luma.gl/core';
import { Deck, DeckGLRef, PickingInfo } from 'deck.gl/typed';
import _ from 'lodash';
import { useCallback, useEffect, useMemo } from 'react';
import { useRecoilCallback, useSetRecoilState } from 'recoil';

import { ViewLayer } from '@/lib/data-map/view-layers';
import { RecoilStateFamily } from '@/lib/recoil/types';

import { allowedGroupLayersState, hoverPositionState, hoverState, selectionState } from './interaction-state';

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

  viewLayer: ViewLayer;
  // logicalLayer: string;

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

function processVectorTarget(info: PickingInfo): VectorTarget {
  const { object } = info;

  return object
    ? {
        feature: object,
      }
    : null;
}

function processTargetByType(type: InteractionStyle, info: PickingInfo) {
  return type === 'raster' ? processRasterTarget(info) : processVectorTarget(info);
}

function processPickedObject(
  info: PickingInfo,
  type: InteractionStyle,
  groupName: string,
  viewLayerLookup: Record<string, ViewLayer>,
  lookupViewForDeck: (deckLayerId: string) => string,
) {
  const deckLayerId = info.layer.id;
  const viewLayerId = lookupViewForDeck(deckLayerId);
  const target = processTargetByType(type, info);

  return (
    target && {
      interactionGroup: groupName,
      interactionStyle: type,
      viewLayer: viewLayerLookup[viewLayerId],
      target,
    }
  );
}

function useSetInteractionGroupState(
  stateFamily: RecoilStateFamily<InteractionTarget<any> | InteractionTarget<any>[], string>,
) {
  return useRecoilCallback(({ set }) => {
    return (groupName: string, value: InteractionTarget<any> | InteractionTarget<any>[]) => {
      set(stateFamily(groupName), value);
    };
  });
}

export function useInteractions(
  viewLayers: ViewLayer[],
  lookupViewForDeck: (deckLayerId: string) => string,
  interactionGroups: InteractionGroupConfig[],
) {
  const setHoverXY = useSetRecoilState(hoverPositionState);

  const setInteractionGroupHover = useSetInteractionGroupState(hoverState);
  const setInteractionGroupSelection = useSetInteractionGroupState(selectionState);

  const interactionGroupLookup = useMemo(() => _.keyBy(interactionGroups, 'id'), [interactionGroups]);

  const primaryGroup = interactionGroups[0].id;
  const primaryGroupPickingRadius = interactionGroupLookup[primaryGroup].pickingRadius;

  const interactiveLayers = useMemo(() => viewLayers.filter((x) => x.interactionGroup), [viewLayers]);
  const viewLayerLookup = useMemo(() => _.keyBy(interactiveLayers, (layer) => layer.id), [interactiveLayers]);
  const activeGroups = useMemo(
    () => _.groupBy(interactiveLayers, (viewLayer) => viewLayer.interactionGroup),
    [interactiveLayers],
  );

  const setAllowedGroupLayers = useSetRecoilState(allowedGroupLayersState);

  useEffect(() => {
    setAllowedGroupLayers(_.mapValues(activeGroups, (viewLayers) => viewLayers.map((viewLayer) => viewLayer.id)));
  }, [activeGroups, setAllowedGroupLayers]);

  const onHover = useCallback(
    (info: any, deck: Deck) => {
      const { x, y } = info;

      for (const [groupName, layers] of Object.entries(activeGroups)) {
        const layerIds = layers.map((layer) => layer.id);
        const interactionGroup = interactionGroupLookup[groupName];
        const { type, pickingRadius: radius, pickMultiple } = interactionGroup;

        const pickingParams = { x, y, layerIds, radius };

        if (pickMultiple) {
          const pickedObjects: PickingInfo[] = deck.pickMultipleObjects(pickingParams);
          const interactionTargets: InteractionTarget<any>[] = pickedObjects
            .map((info) => processPickedObject(info, type, groupName, viewLayerLookup, lookupViewForDeck))
            .filter(Boolean);

          setInteractionGroupHover(groupName, interactionTargets);
        } else {
          const info: PickingInfo = deck.pickObject(pickingParams);
          let interactionTarget: InteractionTarget<any> =
            info && processPickedObject(info, type, groupName, viewLayerLookup, lookupViewForDeck);

          setInteractionGroupHover(groupName, interactionTarget);
        }
      }

      setHoverXY([x, y]);
    },
    [activeGroups, lookupViewForDeck, interactionGroupLookup, setHoverXY, setInteractionGroupHover, viewLayerLookup],
  );

  const onClick = useCallback(
    (info: any, deck: DeckGLRef) => {
      const { x, y } = info;
      for (const [groupName, viewLayers] of Object.entries(activeGroups)) {
        const viewLayerIds = viewLayers.map((layer) => layer.id);

        const interactionGroup = interactionGroupLookup[groupName];
        const { type, pickingRadius: radius } = interactionGroup;

        // currently only supports selecting vector features
        if (interactionGroup.type === 'vector') {
          const info = deck.pickObject({ x, y, layerIds: viewLayerIds, radius });
          let selectionTarget = info && processPickedObject(info, type, groupName, viewLayerLookup, lookupViewForDeck);

          setInteractionGroupSelection(groupName, selectionTarget);
        }
      }
    },
    [activeGroups, lookupViewForDeck, interactionGroupLookup, setInteractionGroupSelection, viewLayerLookup],
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

  const layerFilter = ({ layer: deckLayer, renderPass }) => {
    if (renderPass === 'picking:hover') {
      const viewLayerId = lookupViewForDeck(deckLayer.id);
      const interactionGroup = viewLayerId && viewLayerLookup[viewLayerId]?.interactionGroup;

      return interactionGroup ? hoverPassGroups.has(interactionGroup) : false;
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
