import { useCallback, useState } from 'react';

const networkLayers = [
  'elec_edges_high',
  'elec_edges_low',
  'elec_nodes',
  'rail_edges',
  'rail_nodes',
  'road_edges',
  'road_bridges',
  'water_potable_edges',
  'water_potable_nodes',
];

export const useNetworkSelection = () => {
  const [networkSelection, setNetworkSelection] = useState(Object.fromEntries(networkLayers.map((nl) => [nl, false])));

  const updateNetworkSelection = useCallback(
    (networkSelectionUpdate) => {
      setNetworkSelection({ ...networkSelection, ...networkSelectionUpdate });
    },
    [networkSelection],
  );

  return {
    networkSelection,
    setNetworkSelection: updateNetworkSelection,
    networkVisibilitySet: networkSelection,
  };
};
