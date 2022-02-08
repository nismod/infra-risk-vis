import DeckGL from 'deck.gl';
import { readPixelsToArray } from '@luma.gl/core';
import _ from 'lodash';
import { useCallback, useMemo } from 'react';
import { useRecoilCallback, useSetRecoilState } from 'recoil';

import { VIEW_LAYERS } from 'config/view-layers';

import { hoverState, hoverPositionState } from './interaction-state';

export interface InteractionGroupConfig {
  id: string;
  type: 'vector' | 'raster';
  options?: {
    pickingRadius?: number;
  };
}

export interface InteractionTarget<T> {
  interactionGroup: string;
  interactionStyle: string;

  viewLayer: string;
  logicalLayer: string;

  target: T;
}

function processRasterHover(layerId, info): InteractionTarget<any> {
  const { bitmap, sourceLayer } = info;
  if (bitmap) {
    const pixelColor = readPixelsToArray(sourceLayer.props.image, {
      sourceX: bitmap.pixel[0],
      sourceY: bitmap.pixel[1],
      sourceWidth: 1,
      sourceHeight: 1,
      sourceType: undefined,
    });

    return pixelColor[3] ? pixelColor : null;
  }
}

const pickingRadius = 8;

export function useInteractions(
  viewLayersSpec: Record<string, any>,
  interactionGroups: Record<string, InteractionGroupConfig>,
) {
  /**
   * {
   *  assets: [elec_edges_high, elec_edges_low, water_irrigation_nodes, water_irrigation_edges],
   *  hazard: [hazard__rp_m, ],
   *  region: ['boundaries-parish', 'boundaries-community', 'boundaries-enumeration']
   * }
   */

  const setHoverXY = useSetRecoilState(hoverPositionState);

  const setInteractionGroupHover = useRecoilCallback(({ set }) => {
    return (interactionGroup: string, target: InteractionTarget<any> | InteractionTarget<any>[]) => {
      set(hoverState(interactionGroup), target);
    };
  });

  const activeGroups = useMemo(
    () => _.groupBy(_.keys(viewLayersSpec), (layerName) => VIEW_LAYERS[layerName].interactionGroup),
    [viewLayersSpec],
  );

  console.log(activeGroups);

  const onHover = useCallback(
    (info: any, deck: DeckGL) => {
      const { x, y } = info;

      for (const [groupName, layerIds] of Object.entries(activeGroups)) {
        if (layerIds.length === 0) {
          setInteractionGroupHover(groupName, null);
          continue;
        }

        const interactionGroup = interactionGroups[groupName];

        if (interactionGroup.type === 'vector') {
          const pickedVectorInfo = deck.pickObject({ x, y, layerIds, radius: interactionGroup.options?.pickingRadius });

          if (pickedVectorInfo) {
            const layerId = pickedVectorInfo.layer.id;
            const { object } = pickedVectorInfo;

            setInteractionGroupHover(groupName, {
              interactionGroup: groupName,
              interactionStyle: interactionGroup.type,
              viewLayer: layerId,
              logicalLayer:
                VIEW_LAYERS[layerId].getLogicalLayer?.({ deckLayerId: layerId, feature: object }) ?? layerId,
              target: object,
            });
          } else {
            setInteractionGroupHover(groupName, null);
          }
        } else if (interactionGroup.type === 'raster') {
          const pickedObjects = deck.pickMultipleObjects({ x, y, layerIds });
          const interactionTargets: InteractionTarget<any>[] = [];
          for (const info of pickedObjects) {
            const layerId = info.layer.id;

            const rasterTarget = processRasterHover(layerId, info);

            if (rasterTarget) {
              interactionTargets.push({
                interactionGroup: groupName,
                interactionStyle: interactionGroup.type,
                viewLayer: layerId,
                logicalLayer:
                  VIEW_LAYERS[layerId].getLogicalLayer?.({ deckLayerId: layerId, feature: null }) ?? layerId,
                target: rasterTarget,
              });
            }
          }
          setInteractionGroupHover(groupName, interactionTargets);
        }
      }

      setHoverXY([x, y]);
    },
    [activeGroups, interactionGroups, setHoverXY, setInteractionGroupHover],
  );

  // const onClick = useCallback(
  //   (info: any, deck: DeckGL) => {
  //     const { x, y } = info;
  //     const pickedVectorInfo = deck.pickObject({ x, y, layerIds: activeGroups, radius: pickingRadius });

  //     if (pickedVectorInfo) {
  //       setSelectedFeature(processVectorHover(pickedVectorInfo.layer.id, pickedVectorInfo));
  //     } else {
  //       setSelectedFeature(null);
  //     }
  //   },
  //   [vectorLayerIds],
  // );

  const layerFilter = ({ layer, renderPass }) => {
    if (renderPass === 'picking:hover') {
      // don't render raster and region layers on hover picking pass (but render them for manual picking)
      if (layer.id.match(/^(coastal|fluvial|surface|cyclone|boundaries)/)) return false;
    }
    return true;
  };

  return {
    onHover,
    // onClick,
    layerFilter,
    pickingRadius: 8,
  };
}
