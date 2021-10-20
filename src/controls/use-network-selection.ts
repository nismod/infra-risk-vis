import { useCallback, useState } from 'react';

const networkLayers = [
  'elec_edges_high',
  'elec_edges_low',
  'elec_nodes_source',
  'elec_nodes_sink',
  'elec_nodes_junction',
  'rail_edges',
  'rail_nodes',
  'road_edges',
  'road_bridges',
  'port_areas',
  'airport_areas',
  'water_potable_edges',
  'water_potable_nodes',
  'water_irrigation_edges',
  'water_irrigation_nodes',
  'water_waste_edges',
  'water_waste_nodes',
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
