import { useCallback, useEffect, useState } from 'react';

import { hazardTypes } from '../../config/data/hazards';

export interface SingleHazardSelection {
  show: boolean;
  paramSelection: {
    returnPeriod?: number;
    epoch?: number;
    rcp?: string;
    confidence?: string | number;
  };
  paramOptions: {
    returnPeriod?: number[];
    epoch?: number[];
    rcp?: string[];
    confidence?: (string | number)[];
  };
}

export type HazardSelectionSet = { [k: string]: SingleHazardSelection };

function emptySelection() {
  return Object.fromEntries(hazardTypes.map((ht) => [ht, false]));
}

export const useHazardSelection = (forceSingle = false) => {
  const [hazardSelection, setHazardSelection] = useState(emptySelection());

  useEffect(() => {
    if (forceSingle) {
      let selectedKeys = Object.entries(hazardSelection)
        .filter(([key, value]) => value)
        .map(([key]) => key);
      if (selectedKeys.length > 1) {
        const firstKey = selectedKeys[0];
        setHazardSelection({ ...emptySelection(), [firstKey]: true });
      }
    }
  }, [forceSingle, hazardSelection]);

  const updateHazardSelection = useCallback(
    (hazardType: string, show: boolean) => {
      const base = forceSingle ? emptySelection() : hazardSelection;
      setHazardSelection({ ...base, [hazardType]: show });
    },
    [forceSingle, hazardSelection],
  );

  const clearSelection = useCallback(() => {
    setHazardSelection(emptySelection());
  }, []);

  return {
    hazardSelection,
    setSingleHazardShow: updateHazardSelection,
    clearSelection,
  };
};
