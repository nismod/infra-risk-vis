import { useCallback, useState } from 'react';

import { hazardConfig, hazardTypes } from '../../config/data/hazards';

export const useHazardParams = () => {
  const [hazardParams, setHazardParams] = useState(
    Object.fromEntries(hazardTypes.map((ht) => [ht, hazardConfig[ht].paramDefaults])),
  );

  const [hazardOptions, setHazardOptions] = useState(
    Object.fromEntries(hazardTypes.map((ht) => [ht, hazardConfig[ht].paramDomains])),
  );

  const updatedHazardParam = useCallback(
    (hazardType: string, paramName: string, paramValue: any) => {
      const { paramDomains, paramDependencies = {} } = hazardConfig[hazardType];
      const oldParams = hazardParams[hazardType];

      const newSingleHazardParams = { ...oldParams, [paramName]: paramValue };
      const newSingleHazardOptions = {};

      for (const [param, paramValue] of Object.entries(newSingleHazardParams)) {
        const newParamOptions = paramDependencies[param]?.(newSingleHazardParams) ?? paramDomains[param];

        // if the new options don't include the current param value, switch value to the first option
        if (!newParamOptions.includes(paramValue)) {
          newSingleHazardParams[param] = newParamOptions[0];
        }

        newSingleHazardOptions[param] = newParamOptions;
      }

      setHazardParams({ ...hazardParams, [hazardType]: newSingleHazardParams });
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
