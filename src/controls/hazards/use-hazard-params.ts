import _ from 'lodash';
import { useCallback, useMemo } from 'react';

import { hazardConfig } from '../../config/data/hazards';
import { useDataParams } from '../use-data-params';

export const useHazardParams = () => {
  const fluvialState = useDataParams(hazardConfig.fluvial);
  const surfaceState = useDataParams(hazardConfig.surface);
  const coastalState = useDataParams(hazardConfig.coastal);
  const cycloneState = useDataParams(hazardConfig.cyclone);

  const state = useMemo(
    () => ({
      fluvial: fluvialState,
      surface: surfaceState,
      coastal: coastalState,
      cyclone: cycloneState,
    }),
    [fluvialState, surfaceState, coastalState, cycloneState],
  );

  const updateParam = useCallback(
    (hazardType: string, paramName: string, paramValue: any) => {
      state[hazardType].updateParam(paramName, paramValue);
    },
    [state],
  );

  const hazardParams = _.mapValues(state, (x) => x.params);
  const hazardOptions = _.mapValues(state, (x) => x.options);

  return {
    hazardParams,
    hazardOptions,
    setSingleHazardParam: updateParam,
  };
};
