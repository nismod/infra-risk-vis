import { useMemo } from 'react';

import { hazardTypes } from '../../config/data/hazards';
import { getHazardId } from '../../config/layers';

export const useHazardVisibility = (hazardSelection, hazardParams) => {
  return useMemo(() => {
    const visibility: any = {};

    for (const hazardType of hazardTypes) {
      if (hazardSelection[hazardType]) {
        visibility[getHazardId({ ...hazardParams[hazardType], hazardType })] = true;
      }
    }
    return visibility;
  }, [hazardSelection, hazardParams]);
};
