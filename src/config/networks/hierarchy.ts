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
        label: 'Transmission',
        children: [
          {
            id: 'elec_edges_high',
            label: 'High Voltage',
          },
          {
            id: 'elec_edges_low',
            label: 'Low Voltage',
          },
          {
            id: 'elec_nodes_pole',
            label: 'Poles',
          },
          {
            id: 'elec_nodes_substation',
            label: 'Substations',
          },
        ],
      },
      {
        id: 'elec_nodes_source',
        label: 'Generation',
        children: [
          {
            id: 'elec_nodes_diesel',
            label: 'Diesel',
          },
          {
            id: 'elec_nodes_gas',
            label: 'Gas',
          },
          {
            id: 'elec_nodes_hydro',
            label: 'Hydro',
          },
          {
            id: 'elec_nodes_solar',
            label: 'Solar',
          },
          {
            id: 'elec_nodes_wind',
            label: 'Wind',
          },
        ],
      },
      {
        id: 'elec_nodes_demand',
        label: 'Demand',
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
            id: 'rail_stations',
            label: 'Stations',
          },
          {
            id: 'rail_junctions',
            label: 'Junctions',
          },
        ],
      },
      {
        id: 'road-network',
        label: 'Road Network',
        children: [
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
            id: 'road_bridges',
            label: 'Bridges',
          },
        ],
      },
      {
        id: 'port_areas',
        label: 'Ports',
        children: [
          {
            id: 'port_areas_break',
            label: 'Break',
          },
          {
            id: 'port_areas_container',
            label: 'Container',
          },
          {
            id: 'port_areas_industry',
            label: 'Industry',
          },
          {
            id: 'port_areas_silo',
            label: 'Silo',
          },
        ],
      },
      {
        id: 'air',
        label: 'Airports',
        children: [
          {
            id: 'airport_runways',
            label: 'Runways',
          },
          {
            id: 'airport_terminals',
            label: 'Terminals',
          },
        ],
      },
    ],
  },
];
