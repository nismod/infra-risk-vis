import _ from 'lodash';
import { atomFamily, useRecoilTransaction_UNSTABLE } from 'recoil';

import { DataParamGroupConfig, Param, ParamDomain, resolveParamDependencies } from '@/lib/controls/data-params';
import { toDictionary } from '@/lib/helpers';
import { groupedFamily } from '@/lib/recoil/grouped-family';

import { adaptationDomainsConfig } from '@/config/domains/adaptation';
import { totalDamagesConfig } from '@/config/domains/total-damages';
import { HAZARD_DOMAINS } from '@/config/hazards/domains';

export type DataParamParam = Readonly<{
  group: string;
  param: string;
}>;

export const dataParamConfig: Record<string, DataParamGroupConfig> = {
  ...HAZARD_DOMAINS,
  all: totalDamagesConfig,
  adaptation: adaptationDomainsConfig,
};

export const dataParamNamesByGroup = _.mapValues(dataParamConfig, (groupConfig) => _.keys(groupConfig.paramDefaults));

const dataParamDefaultsByGroup = _.mapValues(dataParamConfig, (groupConfig) =>
  resolveParamDependencies(groupConfig.paramDefaults, groupConfig),
);

export const dataParamState = atomFamily<Param, DataParamParam>({
  key: 'dataParamState',
  default: ({ group, param }: DataParamParam) => dataParamDefaultsByGroup[group][0][param],
});

export const dataParamOptionsState = atomFamily<ParamDomain, DataParamParam>({
  key: 'dataParamOptionsState',
  default: ({ group, param }: DataParamParam) => dataParamDefaultsByGroup[group][1][param],
});

const dataParamNamesByGroupState = atomFamily({
  key: 'dataParamNamesByGroupState',
  default: (group: string) => dataParamNamesByGroup[group],
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
        const paramNames = dataParamNamesByGroup[group];
        const groupParams = toDictionary(
          paramNames,
          (param) => param,
          (param) => get(dataParamState({ group, param })),
        );
        const groupConfig = dataParamConfig[group];

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
