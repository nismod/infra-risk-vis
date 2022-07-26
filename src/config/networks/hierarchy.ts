import { TreeNode } from 'lib/controls/checkbox-tree/tree-node';

interface NetworkLayerData {
  label: string;
}

export const NETWORK_LAYERS_HIERARCHY: TreeNode<NetworkLayerData>[] = [
  {
    id: 'power',
    label: 'Power',
    children: [
      {
        id: 'elec_nodes_source',
        label: 'Generation',
      },
      {
        id: 'elec_edges',
        label: 'Transmission Lines',
      },
      {
        id: 'elec_areas_demand',
        label: 'Demand',
      },
    ],
  },
  {
    id: 'transport',
    label: 'Road Transport',
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
];
