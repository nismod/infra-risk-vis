import { useCallback, useState } from 'react';

type Params = object;
type ParamDomain<PT> = PT[];

export type ParamDomains<PT extends Params> = {
  [K in keyof PT]: ParamDomain<PT[K]>;
};

export type ParamDependencies<PT extends Params> = {
  [K in keyof PT]?: (params: PT) => ParamDomain<PT[K]>;
};

export interface DataParamConfig<P extends Params> {
  paramDomains: ParamDomains<P>;
  paramDefaults: P;
  paramDependencies?: ParamDependencies<P>;
}

function getNewParamOptions(updatedParams, domain, dependencyFn) {
  return dependencyFn?.(updatedParams) ?? domain;
}

function resolveDependencies<P extends object>(
  updatedParams: P,
  paramDomains: ParamDomains<P>,
  paramDependencies: ParamDependencies<P>,
): [P, ParamDomains<P>] {
  const resolvedParams = { ...updatedParams };
  const newOptions: ParamDomains<P> = {} as any;

  for (const [param, paramValue] of Object.entries(updatedParams)) {
    const newParamOptions = getNewParamOptions(updatedParams, paramDomains[param], paramDependencies[param]);

    // if the new options don't include the current param value, switch value to the first option
    if (!newParamOptions.includes(paramValue)) {
      resolvedParams[param] = newParamOptions[0];
    }

    newOptions[param] = newParamOptions;
  }
  return [resolvedParams, newOptions];
}

export function useDataParams<P extends object>(config: DataParamConfig<P>) {
  const [params, setParams] = useState(config.paramDefaults);
  const [options, setOptions] = useState(config.paramDomains);

  const updateParam = useCallback(
    (paramName: keyof P, paramValue: any) => {
      const { paramDomains, paramDependencies = {} } = config;

      const updatedParams = { ...params, [paramName]: paramValue } as P;
      const [resolvedParams, newOptions] = resolveDependencies(updatedParams, paramDomains, paramDependencies);

      setParams(resolvedParams);
      setOptions(newOptions);
    },
    [config, params],
  );

  return {
    params,
    options,
    updateParam,
  };
}
