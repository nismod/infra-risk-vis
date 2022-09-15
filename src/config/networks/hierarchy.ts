import { TreeNode } from 'lib/controls/checkbox-tree/tree-node';

interface NetworkLayerData {
  label: string;
}

export const NETWORK_LAYERS_HIERARCHY: TreeNode<NetworkLayerData>[] = [
  {
    id: 'rail-network',
    label: 'Rail',
    children: [
      {
        id: 'rail_edges_open',
        label: 'Rail (open)',
      },
      {
        id: 'rail_edges_disused',
        label: 'Rail (disused)',
      },
      {
        id: 'rail_edges_rehabilitation',
        label: 'Rail (rehabilitation)',
      },
      {
        id: 'rail_edges_construction',
        label: 'Rail (construction)',
      },
      {
        id: 'rail_edges_abandoned',
        label: 'Rail (abandoned)',
      },
      {
        id: 'rail_edges_proposed',
        label: 'Rail (proposed)',
      },
      {
        id: 'rail_nodes_station',
        label: 'Station',
      },
      {
        id: 'rail_nodes_stop',
        label: 'Stop',
      },
      {
        id: 'rail_nodes_halt',
        label: 'Halt',
      },
    ],
  },
  {
    id: 'roads',
    label: 'Roads',
    children: [
      {
        id: 'road_edges_motorway',
        label: 'Motorway',
      },
      {
        id: 'road_edges_trunk',
        label: 'Trunk',
      },
      {
        id: 'road_edges_primary',
        label: 'Primary',
      },
      {
        id: 'road_edges_secondary',
        label: 'Secondary',
      },
      {
        id: 'road_edges_tertiary',
        label: 'Tertiary',
      },
    ],
  },
  {
    id: 'port_areas',
    label: 'Ports',
    children: [
      {
        id: 'port_nodes_maritime',
        label: 'Maritime',
      },
      {
        id: 'port_nodes_lake',
        label: 'Lake',
      },
    ]
  },
  {
    id: 'air_nodes',
    label: 'Airports',
  },
];
