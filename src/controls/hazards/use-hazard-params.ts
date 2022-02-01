import _ from 'lodash';
import { useCallback, useState } from 'react';

import { hazardConfig, HazardParams } from '../../config/data/hazards';
import { ParamDependencies, ParamDomains } from '../../config/data/types';

function resolveDependencies<P extends object = HazardParams>(
  updatedParams: P,
  paramDomains: ParamDomains<P>,
  paramDependencies: ParamDependencies<P>,
): [P, ParamDomains<P>] {
  const resolvedParams = { ...updatedParams };
  const newSingleHazardOptions: ParamDomains<P> = {} as any;

  for (const [param, paramValue] of Object.entries(updatedParams)) {
    const newParamOptions = paramDependencies[param]?.(updatedParams) ?? paramDomains[param];

    // if the new options don't include the current param value, switch value to the first option
    if (!newParamOptions.includes(paramValue)) {
      resolvedParams[param] = newParamOptions[0];
    }

    newSingleHazardOptions[param] = newParamOptions;
  }
  return [resolvedParams, newSingleHazardOptions];
}

export const useHazardParams = () => {
  const [hazardParams, setHazardParams] = useState(_.mapValues(hazardConfig, 'paramDefaults'));
  const [hazardOptions, setHazardOptions] = useState(_.mapValues(hazardConfig, 'paramDomains'));

  const updatedHazardParam = useCallback(
    (hazardType: string, paramName: string, paramValue: any) => {
      const { paramDomains, paramDependencies = {} } = hazardConfig[hazardType];
      const oldParams = hazardParams[hazardType];

      const updatedSingleHazardParams = { ...oldParams, [paramName]: paramValue };
      const [resolvedSingleHazardParams, newSingleHazardOptions] = resolveDependencies(
        updatedSingleHazardParams,
        paramDomains,
        paramDependencies,
      );

      setHazardParams({ ...hazardParams, [hazardType]: resolvedSingleHazardParams });
      setHazardOptions({ ...hazardOptions, [hazardType]: newSingleHazardOptions });
    },
    [hazardOptions, hazardParams],
  );

  return {
    hazardParams,
    hazardOptions,
    setSingleHazardParam: updatedHazardParam,
  };
};
