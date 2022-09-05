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
        id: 'rail_edges',
        label: 'Railways',
      },
      {
        id: 'rail_stations',
        label: 'Stations',
      },
    ],
  },
  {
    id: 'roads',
    label: 'Roads',
    children: [
      {
        id: 'road_edges_class_a',
        label: 'Class A',
      },
      {
        id: 'road_edges_class_b',
        label: 'Class B',
      },
      {
        id: 'road_edges_class_c',
        label: 'Class C',
      },
      {
        id: 'road_edges_metro',
        label: 'Metro',
      },
      {
        id: 'road_edges_track',
        label: 'Track',
      },
      {
        id: 'road_edges_other',
        label: 'Other',
      },
    ],
  },
  {
    id: 'port_areas',
    label: 'Ports',
  },
  {
    id: 'air',
    label: 'Airports',
  },
];
