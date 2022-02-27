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
        id: 'power-lines',
        label: 'Power Lines',
        children: [
          {
            id: 'elec_edges_high',
            label: 'High Voltage',
          },
          {
            id: 'elec_edges_low',
            label: 'Low Voltage',
          },
        ],
      },
      {
        id: 'power-nodes',
        label: 'Power Nodes',
        children: [
          {
            id: 'elec_nodes_source',
            label: 'Generation',
          },
          {
            id: 'elec_nodes_sink',
            label: 'Demand',
          },
          {
            id: 'elec_nodes_junction',
            label: 'Junctions',
          },
        ],
      },
    ],
  },
  {
    id: 'transport',
    label: 'Transport',
    children: [
      {
        id: 'rail-network',
        label: 'Rail Network',
        children: [
          {
            id: 'rail_edges',
            label: 'Railways',
          },
          {
            id: 'rail_nodes',
            label: 'Stations',
          },
        ],
      },
      {
        id: 'road-network',
        label: 'Road Network',
        children: [
          {
            id: 'road_edges',
            label: 'Roads',
          },
          {
            id: 'road_bridges',
            label: 'Bridges',
          },
        ],
      },
      {
        id: 'port_areas',
        label: 'Ports',
      },
      {
        id: 'airport_areas',
        label: 'Airports',
      },
    ],
  },
  {
    id: 'water',
    label: 'Water',
    children: [
      {
        id: 'water-supply',
        label: 'Water Supply Network',
        children: [
          {
            id: 'water_potable_edges',
            label: 'Supply Pipelines',
          },
          {
            id: 'water_potable_nodes',
            label: 'Supply Facilities',
          },
        ],
      },
      {
        id: 'water-irrigation',
        label: 'Irrigation Network',
        children: [
          {
            id: 'water_irrigation_edges',
            label: 'Irrigation Canals',
          },
          {
            id: 'water_irrigation_nodes',
            label: 'Irrigation Facilities',
          },
        ],
      },
      {
        id: 'water-waste',
        label: 'Wastewater Network',
        children: [
          {
            id: 'water_waste_edges',
            label: 'Wastewater Pipelines',
          },
          {
            id: 'water_waste_nodes',
            label: 'Wastewater Facilities',
          },
        ],
      },
    ],
  },
  {
    id: 'buildings',
    label: 'Buildings',
  },
];
