import { HAZARD_DOMAINS } from 'config/hazards/domains';
import { totalDamagesConfig } from 'config/data/total-damages';
import { DataParamGroupConfig, Param, ParamDomain, resolveParamDependencies } from 'lib/controls/data-params';
import { toDictionary } from 'lib/helpers';
import { groupedFamily } from 'lib/recoil/grouped-family';
import _ from 'lodash';
import { atomFamily, selectorFamily, useRecoilTransaction_UNSTABLE } from 'recoil';

export type DataParamParam = Readonly<{
  group: string;
  param: string;
}>;

const dataParamConfig: Record<string, DataParamGroupConfig> = {
  ...HAZARD_DOMAINS,
  'total-damages': totalDamagesConfig,
};

export const dataParamGroupConfigState = atomFamily<DataParamGroupConfig, string>({
  key: 'dataParamGroupConfigState',
  default: (group) => dataParamConfig[group],
});

export const dataParamNamesByGroupState = selectorFamily<string[], string>({
  key: 'dataParamNamesByGroupState',
  get:
    (group) =>
    ({ get }) =>
      Object.keys(get(dataParamGroupConfigState(group)).paramDefaults),
});

// TODO: remove duplication of default calculation between dataParamState and dataParamOptionsState
export const dataParamState = atomFamily<Param, DataParamParam>({
  key: 'dataParamState',
  default: ({ group, param }) => {
    const groupConfig = dataParamConfig[group];
    const [resolvedParams] = resolveParamDependencies(groupConfig.paramDefaults, groupConfig);
    return resolvedParams[param];
  },
});

export const dataParamOptionsState = atomFamily<ParamDomain, DataParamParam>({
  key: 'dataParamOptionsState',
  default: ({ group, param }) => {
    const groupConfig = dataParamConfig[group];
    const [, resolvedOptions] = resolveParamDependencies(groupConfig.paramDefaults, groupConfig);
    return resolvedOptions[param];
  },
});

export const dataParamsByGroupState = groupedFamily<Param, DataParamParam>(
  'dataParamsByGroupState',
  dataParamState,
  dataParamNamesByGroupState,
  (group, param) => ({ group, param }),
);

export const dataParamOptionsByGroupState = groupedFamily<ParamDomain, DataParamParam>(
  'dataParamOptionsByGroupState',
  dataParamOptionsState,
  dataParamNamesByGroupState,
  (group, param) => ({ group, param }),
);

export function useUpdateDataParam(group: string, paramId: string) {
  return useRecoilTransaction_UNSTABLE(
    ({ get, set }) =>
      (newValue) => {
        // need to construct group params dict manually as useRecoilTransaction doesn't support reading selectors
        const paramNames = _.keys(dataParamConfig[group].paramDefaults);
        const groupParams = toDictionary(
          paramNames,
          (param) => param,
          (param) => get(dataParamState({ group, param })),
        );
        const groupConfig = get(dataParamGroupConfigState(group));

        const [resolvedParams, resolvedOptions] = resolveParamDependencies<Record<string, any>>(
          { ...groupParams, [paramId]: newValue },
          groupConfig,
        );

        _.forEach(resolvedParams, (resolvedParamValue, paramId) => {
          const recoilParam = { group, param: paramId };
          set(dataParamState(recoilParam), resolvedParamValue);
          set(dataParamOptionsState(recoilParam), resolvedOptions[paramId]);
        });
      },
    [group, paramId],
  );
}
