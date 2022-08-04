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
        children: [
          {
            id: 'elec_nodes_biomass',
            label: 'Biomass',
          },
          {
            id: 'elec_nodes_coal',
            label: 'Coal',
          },
          {
            id: 'elec_nodes_gas',
            label: 'Gas',
          },
          {
            id: 'elec_nodes_geothermal',
            label: 'Geothermal',
          },
          {
            id: 'elec_nodes_hydro',
            label: 'Hydro',
          },
          {
            id: 'elec_nodes_oil',
            label: 'Oil',
          },
          {
            id: 'elec_nodes_solar',
            label: 'Solar',
          },
          {
            id: 'elec_nodes_wind',
            label: 'Wind',
          },
        ]
      },
      {
        id: 'elec_edges',
        label: 'Transmission Lines',
      },
    ],
  },
  {
    id: 'transport',
    label: 'Road Transport',
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
];
