import _ from 'lodash';
import { useEffect } from 'react';
import {
  RecoilValue,
  atomFamily,
  selectorFamily,
  useRecoilTransaction_UNSTABLE,
  useRecoilValue,
  useRecoilValueLoadable,
  useSetRecoilState,
} from 'recoil';

import {
  DataParamGroupConfig,
  ParamDomain,
  ParamValue,
  resolveParamDependencies,
} from '@/lib/controls/data-params';

/**
 * Data params config state (domains, defaults, dependencies) per group
 */
export const paramsConfigState = atomFamily<DataParamGroupConfig, string>({
  key: 'paramsConfigState',
  default: () => new Promise(() => {}),
});

interface ValueAndOptions<T = any> {
  value: T;
  options: T[];
}

/**
 * Data params state (current value + available options) per group
 */
export const paramsState = atomFamily<Record<string, ValueAndOptions>, string>({
  key: 'paramsState',
  default: () => new Promise(() => {}),
});

/**
 * Initialize data params state (value/options) and config for a group, based on an external config.
 *
 * `useUpdateDataParams` relies on `useRecoilTransaction` which currently doesn't support reading from selectors.
 * This forces `paramsConfigState` to be an atom family, but that prevents loading the config from the API with async selectors.
 * `useLoadParamsConfig` can be used in a component to load the external config (with Suspense while the config's loading),
 * and only set `paramsState` and `paramsConfigState` when the value is loaded.
 *
 * @param configState a recoil state representing the loaded configuration. Can be async.
 * @param targetGroup string representing the target data params group
 */
export function useLoadParamsConfig(
  configState: RecoilValue<DataParamGroupConfig>,
  targetGroup: string,
) {
  const config = useRecoilValue(configState);
  const setTargetConfig = useSetRecoilState(paramsConfigState(targetGroup));
  const setTargetState = useSetRecoilState(paramsState(targetGroup));
  const { state: loadableState } = useRecoilValueLoadable(paramsConfigState(targetGroup));

  useEffect(() => {
    // only load config once
    if (loadableState !== 'hasValue') {
      const [values, options] = resolveParamDependencies(config.paramDefaults, config);
      const initialState = _.mapValues(values, (value, key) => ({ value, options: options[key] }));
      setTargetState(initialState);
      setTargetConfig(config);
    }
  }, [config, setTargetConfig, setTargetState, loadableState]);
}

/**
 * Returns a recoil transaction callback for updating a data param with a new value.
 *
 * Resolves dependencies between individual parameters, to update the available options for each param.
 * @param group data param group name
 * @param paramId data param name
 * @returns the value update callback
 */
export function useUpdateDataParam(group: string, paramId: string) {
  return useRecoilTransaction_UNSTABLE(
    ({ get, set }) =>
      (newValue) => {
        const config = get(paramsConfigState(group));

        const state = get(paramsState(group));

        const oldValues = _.mapValues(state, (x) => x.value);
        const newValues = { ...oldValues, [paramId]: newValue };

        const [resolvedValues, resolvedOptions] = resolveParamDependencies<Record<string, any>>(
          newValues,
          config,
        );

        const resolvedState = _.mapValues(state, (x, key) => ({
          value: resolvedValues[key],
          options: resolvedOptions[key],
        }));

        set(paramsState(group), resolvedState);
      },
    [group, paramId],
  );
}

/**
 * Recoil atom/selector family param describing a data parameter (data group + data param name)
 */
export type DataParamParam = Readonly<{
  group: string;
  param: string;
}>;

/**
 * Data param current value per group+param
 */
export const paramValueState = selectorFamily<ParamValue, DataParamParam>({
  key: 'paramValueState',
  get:
    ({ group, param }) =>
    ({ get }) =>
      get(paramsState(group))?.[param]?.value,
});

/**
 * Data param available options per group+param
 */
export const paramOptionsState = selectorFamily<ParamDomain, DataParamParam>({
  key: 'paramOptionsState',
  get:
    ({ group, param }) =>
    ({ get }) =>
      get(paramsState(group))?.[param]?.options ?? [],
});

/**
 * Dictionary of all param current values, per group
 */
export const dataParamsByGroupState = selectorFamily({
  key: 'dataParamsByGroupState',
  get:
    (group: string) =>
    ({ get }) =>
      _.mapValues(get(paramsState(group)), (x) => x.value),
});
