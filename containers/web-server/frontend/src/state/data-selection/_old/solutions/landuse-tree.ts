import _ from 'lodash';
import { atom, selector } from 'recoil';

import { CheckboxTreeState, buildTreeConfig } from '@/lib/controls/checkbox-tree/CheckboxTree';

import { LandUseOption } from '@/config/_old/solutions/domains';
import { LANDUSE_HIERARCHY } from '@/config/_old/solutions/landuse-hierarchy';

export const landuseTreeExpandedState = atom<string[]>({
  key: 'landuseTreeExpandedState',
  default: [],
});

export const landuseTreeConfig = buildTreeConfig(LANDUSE_HIERARCHY);

export const landuseTreeCheckboxState = atom<CheckboxTreeState>({
  key: 'landuseTreeCheckboxState',
  default: {
    checked: _.mapValues(landuseTreeConfig.nodes, () => true),
    indeterminate: _.mapValues(landuseTreeConfig.nodes, () => false),
  },
});

export const landuseFilterState = selector<Record<LandUseOption, boolean>>({
  key: 'landuseFilterState',
  get: ({ get }) => {
    const checkboxState = get(landuseTreeCheckboxState).checked;

    return _.pickBy(checkboxState, (value, id) => !landuseTreeConfig.nodes[id].children) as Record<
      LandUseOption,
      boolean
    >;
  },
});
