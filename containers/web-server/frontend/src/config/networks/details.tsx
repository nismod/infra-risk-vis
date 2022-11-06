import { List } from '@mui/material';
import { FC } from 'react';

import { numFormat, numRangeFormat, paren } from '@/lib/helpers';
import { DataItem } from '@/lib/ui/data-display/DataItem';

import {
  DetailHeader,
  DetailsComponentProps,
  DetailsComponentType,
  IdSubheader,
} from '@/details/features/detail-components';

import { NetworkLayerType } from './metadata';

export const AirportDetails: FC<DetailsComponentProps> = ({ f }) => (
  <>
    <DetailHeader>{f.name}</DetailHeader>
    <IdSubheader id={f.asset_id} />
    <List>
      <DataItem label="Passengers (people/year)" value={f.passenger_number} />
      <DataItem label="Freight (tonnes/year)" value={f.freight_tonnes} />
      <DataItem
        label={`Rehabilitation cost (${f.cost_unit})`}
        value={`${numFormat(f.cost_mean)} ${paren(numRangeFormat(f.cost_min, f.cost_max))}`}
      />
    </List>
  </>
);

export const PowerLineDetails: FC<DetailsComponentProps> = ({ f }) => (
  <>
    {f.name && <DetailHeader>{f.name}</DetailHeader>}
    <IdSubheader id={f.asset_id} />
    <List>
      <DataItem label="Connection" value={`${f.from_type} ${f.from_id}–${f.to_type} ${f.to_id} `} />
      <DataItem label="Voltage (kV)" value={f.voltage_kV} />
      <DataItem label="Length (m)" value={f.length} />
      <DataItem
        label={`Rehabilitation cost (${f.cost_unit})`}
        value={`${numFormat(f.cost_mean)} (${numFormat(f.cost_min)}–${numFormat(f.cost_max)})`}
      />
    </List>
  </>
);

export const PowerGenerationNodeDetails: FC<DetailsComponentProps> = ({ f }) => (
  <>
    <DetailHeader>{f.title}</DetailHeader>
    <IdSubheader id={f.asset_id} />
    <List>
      <DataItem label="Capacity (MW)" value={f.capacity} />
      <DataItem label={`Energy intensity (${f.ei_uom})`} value={f.ei} />
      <DataItem
        label={`Rehabilitation cost (${f.cost_unit})`}
        value={`${numFormat(f.cost_mean)} (${numFormat(f.cost_min)}–${numFormat(f.cost_max)})`}
      />
    </List>
  </>
);

export const PowerDemandNodeDetails: FC<DetailsComponentProps> = ({ f }) => (
  <>
    <IdSubheader id={f.asset_id} />
    <List>
      <DataItem label="Population served" value={f.population} />
      <DataItem label={`Energy intensity (${f.ei_uom})`} value={f.ei} />
    </List>
  </>
);

export const PowerJunctionNodeDetails: FC<DetailsComponentProps> = ({ f }) => (
  <>
    <IdSubheader id={f.asset_id} />
    <List>
      <DataItem label={`Energy intensity (${f.ei_uom})`} value={f.ei} />
      <DataItem
        label={`Rehabilitation cost (${f.cost_unit})`}
        value={`${numFormat(f.cost_mean)} (${numFormat(f.cost_min)}–${numFormat(f.cost_max)})`}
      />
    </List>
  </>
);

export const IrrigationDetails: FC<DetailsComponentProps> = ({ f }) => (
  <>
    <IdSubheader id={f.asset_id} />
    <List>
      <DataItem label={`Rehabilitation cost (${f.cost_unit})`} value={numRangeFormat(f.cost_min, f.cost_max)} />
      <DataItem label="Notes" value={f.comment} />
    </List>
  </>
);

export const WaterPipelineDetails: FC<DetailsComponentProps> = ({ f }) => (
  <>
    <IdSubheader id={f.asset_id} />
    <List>
      <DataItem label="Material" value={f.Material} />
      <DataItem label="Diameter (m)" value={f.Diameter} />
      <DataItem label="Length (m)" value={f.Length} />
      <DataItem label={`Rehabilitation cost (${f.cost_unit})`} value={numRangeFormat(f.cost_min, f.cost_max)} />
      <DataItem label="Notes" value={f.comment} />
    </List>
  </>
);

export const PortDetails: FC<DetailsComponentProps> = ({ f }) => (
  <>
    <DetailHeader>{f.name}</DetailHeader>
    <IdSubheader id={f.asset_id} />
    <List>
      <DataItem label="Vessels (vessels/year)" value={f.vessels_number} />
      <DataItem label="Passengers (people/year)" value={f.passenger_number} />
      <DataItem label="Commodity" value={f.commodity} />
      <DataItem label="Export (tonnes/year)" value={f.export_tonnes} />
      <DataItem label="Import (tonnes/year)" value={f.import_tonnes} />
      <DataItem label="Import of vehicles (tonnes/year)" value={f.import_vehicles_tonnes} />
      <DataItem label="Transhipment (tonnes/year)" value={f.transhipment_tonnes} />
      <DataItem
        label={`Rehabilitation cost (${f.cost_unit})`}
        value={`${numFormat(f.cost_mean)} ${paren(numRangeFormat(f.cost_min, f.cost_max))}`}
      />
      <DataItem label="Notes" value={f.comment} />
    </List>
  </>
);

export const WaterSupplyNodeDetails: FC<DetailsComponentProps> = ({ f }) => (
  <>
    <DetailHeader>{f.name}</DetailHeader>
    <IdSubheader id={f.asset_id} />
    <List>
      <DataItem label="Tank" value={f.Tank_Type} />
      <DataItem label="Population served" value={f.asset_pop_new_2010} />
      <DataItem label="Capacity (mgd, millions of gallons/day)" value={f['capacity (mgd)']} />
      <DataItem label={`Rehabilitation cost (${f.cost_unit})`} value={numRangeFormat(f.cost_min, f.cost_max)} />
      <DataItem label="Notes" value={f.comment} />
    </List>
  </>
);

export const RailEdgeDetails: FC<DetailsComponentProps> = ({ f }) => (
  <>
    <DetailHeader>{f.rail_sect}</DetailHeader>
    <IdSubheader id={f.asset_id} />
    <List>
      <DataItem label="Connection" value={`${f.from_node}–${f.to_node}`} />
      <DataItem label="Owner" value={f.type} />
      <DataItem label="Status" value={f.status} />
      <DataItem label="Length (m)" value={f.shape_length} />
      <DataItem
        label={`Rehabilitation cost (${f.cost_unit})`}
        value={`${numFormat(f.cost_mean)} (${numFormat(f.cost_min)}–${numFormat(f.cost_max)})`}
      />
    </List>
  </>
);

export const RailNodeDetails: FC<DetailsComponentProps> = ({ f }) => (
  <>
    <DetailHeader>{f.Station}</DetailHeader>
    <IdSubheader id={f.asset_id} />
    <List>
      <DataItem label="Lines" value={f.Lines} />
      <DataItem label="Status" value={`${f.status} ${paren(f.Condition)}`} />
      <DataItem
        label={`Rehabilitation cost (${f.cost_unit})`}
        value={`${numFormat(f.cost_mean)} (${numFormat(f.cost_min)}–${numFormat(f.cost_max)})`}
      />
    </List>
  </>
);

export const RoadEdgeDetails: FC<DetailsComponentProps> = ({ f }) => (
  <>
    <DetailHeader>
      {f.street_name}
      {f.street_name ? ', ' : ''}
      {f.section_name}
    </DetailHeader>
    <IdSubheader id={f.asset_id} />
    <List>
      <DataItem label="Connection" value={`${f.from_node}–${f.to_node}`} />
      <DataItem label="Street type" value={`${f.street_type ? f.street_type : 'none'}`} />
      <DataItem label="Construction" value={f.road_construction} />
      <DataItem label="Length (m)" value={f.length_m} />
      <DataItem label="Width (m)" value={f.road_width} />
      <DataItem label="Vertical alignment" value={f.vertalignm} />
      <DataItem label="Traffic (vehicles/day)" value={f.traffic_count} />
      <DataItem
        label={`Rehabilitation cost (${f.cost_unit})`}
        value={`${numFormat(f.cost_mean)} (${numFormat(f.cost_min)}–${numFormat(f.cost_max)})`}
      />
      <DataItem label={`Reopening cost (${f.cost_reopen_unit})`} value={`${numFormat(f.cost_reopen)}`} />
    </List>
  </>
);

export const RoadJunctionDetails: FC<DetailsComponentProps> = ({ f }) => (
  <>
    <DetailHeader>{f.name}</DetailHeader>
    <IdSubheader id={f.asset_id} />
    <List>
      <DataItem
        label={`Rehabilitation cost (${f.cost_unit})`}
        value={`${numFormat(f.cost_mean)} ${paren(numRangeFormat(f.cost_min, f.cost_max))}`}
      />
      <DataItem label={`Reopening cost (${f.cost_unit})`} value={`${numFormat(f.cost_reopen)}`} />
    </List>
  </>
);

export const BridgeDetails: FC<DetailsComponentProps> = ({ f }) => (
  <>
    <IdSubheader id={f.asset_id} />
    <List>
      <DataItem label="Source ID" value={f.BRIDGEID} />
      <DataItem
        label={`Rehabilitation cost (${f.cost_unit})`}
        value={`${numFormat(f.cost_mean)} ${paren(numRangeFormat(f.cost_min, f.cost_max))}`}
      />
      <DataItem label={`Reopening cost (${f.cost_unit})`} value={`${numFormat(f.cost_reopen)}`} />
    </List>
  </>
);

export const WastewaterNodeDetails: FC<DetailsComponentProps> = ({ f }) => (
  <>
    <DetailHeader>{f.Name}</DetailHeader>
    <IdSubheader id={f.asset_id} />
    <List>
      <DataItem label="Treatment type" value={f.Type_of_Tr} />
      <DataItem label="Status" value={f.Status} />
      <DataItem label="Capacity (mgd, millions of gallons/day)" value={f['capacity (mgd)']} />
      <DataItem label={`Rehabilitation cost (${f.cost_unit})`} value={numRangeFormat(f.cost_min, f.cost_max)} />
      <DataItem label="Notes" value={f.comment} />
    </List>
  </>
);

export const INFRASTRUCTURE_LAYER_DETAILS: Record<NetworkLayerType, DetailsComponentType> = {
  /*
      airport_terminals: AirportDetails,
      airport_runways: AirportDetails,

      port_areas_break: PortDetails,
      port_areas_container: PortDetails,
      port_areas_industry: PortDetails,
      port_areas_silo: PortDetails,

      */
  rail_edges: RailEdgeDetails,
  // rail_stations: RailNodeDetails,
  rail_nodes: RailNodeDetails,

  // road_bridges: BridgeDetails,
  road_edges_motorway: RoadEdgeDetails,
  road_edges_trunk: RoadEdgeDetails,
  road_edges_primary: RoadEdgeDetails,
  road_edges_secondary: RoadEdgeDetails,
  road_edges_tertiary: RoadEdgeDetails,

  /*
      water_irrigation_edges: IrrigationDetails,
      water_irrigation_nodes: IrrigationDetails,

      water_potable_nodes_booster: WaterSupplyNodeDetails,
      water_potable_nodes_catchment: WaterSupplyNodeDetails,
      water_potable_nodes_entombment: WaterSupplyNodeDetails,
      water_potable_nodes_filter: WaterSupplyNodeDetails,
      water_potable_nodes_intake: WaterSupplyNodeDetails,
      water_potable_nodes_well: WaterSupplyNodeDetails,
      water_potable_nodes_pump: WaterSupplyNodeDetails,
      water_potable_nodes_relift: WaterSupplyNodeDetails,
      water_potable_nodes_reservoir: WaterSupplyNodeDetails,
      water_potable_nodes_river_source: WaterSupplyNodeDetails,
      water_potable_nodes_spring: WaterSupplyNodeDetails,
      water_potable_nodes_tank: WaterSupplyNodeDetails,
      water_potable_nodes_sump: WaterSupplyNodeDetails,
      water_potable_nodes_tp: WaterSupplyNodeDetails,
      water_potable_edges: WaterPipelineDetails,

      water_waste_nodes_sump: WastewaterNodeDetails,
      water_waste_nodes_pump: WastewaterNodeDetails,
      water_waste_nodes_relift: WastewaterNodeDetails,
      water_waste_nodes_wwtp: WastewaterNodeDetails,
      water_waste_sewer_gravity: WaterPipelineDetails,
      water_waste_sewer_pressure: WaterPipelineDetails,
      */

  elec_edges_high: PowerLineDetails,
  elec_edges_low: PowerLineDetails,
  elec_nodes_pole: PowerJunctionNodeDetails,
  elec_nodes_substation: PowerJunctionNodeDetails, // TODO create own component

  elec_nodes_demand: PowerDemandNodeDetails,

  elec_nodes_diesel: PowerGenerationNodeDetails,
  elec_nodes_gas: PowerGenerationNodeDetails,
  elec_nodes_hydro: PowerGenerationNodeDetails,
  elec_nodes_solar: PowerGenerationNodeDetails,
  elec_nodes_wind: PowerGenerationNodeDetails,
};
