import { atomFamily, selector, selectorFamily, useRecoilTransaction_UNSTABLE } from 'recoil';

import { ViewLayer, ViewLayerParams } from 'lib/data-map/view-layers';

import { viewLayersFlatState } from './view-layers-flat';
import _ from 'lodash';
import { selectionState } from 'lib/data-map/interactions/interaction-state';
import { styleParamsState } from 'state/style-params';

export const viewLayerState = atomFamily({
  key: 'viewLayerState',
  default: null,
});

export const useSaveViewLayers = () => {
  return useRecoilTransaction_UNSTABLE(
    ({ set }) =>
      (viewLayers: ViewLayer[]) =>
        viewLayers.forEach((viewLayer) => set(viewLayerState(viewLayer.id), viewLayer)),
  );
};

export const singleViewLayerParamsState = selectorFamily<ViewLayerParams, string>({
  key: 'singleViewLayerParamsState',
  get:
    (viewLayerId: string) =>
    ({ get }) => {
      const viewLayer = get(viewLayerState(viewLayerId));

      if (viewLayer?.group === 'infrastructure') {
        const interactionGroup = viewLayer.interactionGroup;
        const groupSelection = get(selectionState(interactionGroup));

        const styleParams = get(styleParamsState);

        return {
          styleParams,
          selection: groupSelection?.viewLayer.id === viewLayer.id ? groupSelection : null,
        };
      } else return {};
    },
});

export const viewLayersParamsState = selector<Record<string, ViewLayerParams>>({
  key: 'viewLayersParamsState',
  get: ({ get }) => {
    const viewLayers = get(viewLayersFlatState);

    return _(viewLayers)
      .keyBy('id')
      .mapValues((viewLayer) => get(singleViewLayerParamsState(viewLayer.id)))
      .value();
  },
});
