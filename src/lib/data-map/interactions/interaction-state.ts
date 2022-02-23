import _ from 'lodash';
import { atom, atomFamily, selector } from 'recoil';

import { isReset } from 'lib/recoil/is-reset';

import { InteractionTarget } from './use-interactions';

type IT = InteractionTarget<any> | InteractionTarget<any>[];

export function hasHover(target: IT) {
  if (Array.isArray(target)) {
    return target.length > 0;
  }
  return !!target;
}

export const hoverState = atomFamily<IT, string>({
  key: 'hoverState',
  default: null,
});

export const hoverPositionState = atom({
  key: 'hoverPosition',
  default: null,
});

export const selectionState = atomFamily<InteractionTarget<any>, string>({
  key: 'selectionState',
  default: null,
});

type AllowedGroupLayers = Record<string, string[]>;

const allowedGroupLayersImpl = atom<AllowedGroupLayers>({
  key: 'allowedGroupLayersImpl',
  default: {},
});

function filterOneOrArray<T>(items: T | T[], filter: (item: T) => boolean) {
  if (Array.isArray(items)) {
    return items.filter(filter);
  } else {
    return items && filter(items) ? items : null;
  }
}

function filterTargets(oldHoverTargets: IT, allowedLayers: string[]): IT {
  const newLayerFilter = new Set(allowedLayers);
  return filterOneOrArray(oldHoverTargets, (target) => newLayerFilter.has(target.viewLayer.id));
}

export const allowedGroupLayersState = selector<AllowedGroupLayers>({
  key: 'allowedGroupLayersState',
  get: ({ get }) => get(allowedGroupLayersImpl),
  set: ({ get, set, reset }, newAllowedGroups) => {
    const oldAllowedGroupLayers = get(allowedGroupLayersImpl);
    if (isReset(newAllowedGroups)) {
      _.forEach(oldAllowedGroupLayers, (layers, group) => {
        reset(hoverState(group));
        reset(selectionState(group));
      });
    } else {
      for (const group of Object.keys(oldAllowedGroupLayers)) {
        const newAllowedLayers = newAllowedGroups[group];

        if (newAllowedLayers == null || newAllowedLayers.length === 0) {
          reset(hoverState(group));
          reset(selectionState(group));
        } else {
          const oldHoverTargets = get(hoverState(group));
          set(hoverState(group), filterTargets(oldHoverTargets, newAllowedLayers));

          const oldSelectionTargets = get(selectionState(group));
          set(selectionState(group), filterTargets(oldSelectionTargets, newAllowedLayers));
        }
      }
    }

    set(allowedGroupLayersImpl, newAllowedGroups);
  },
});
