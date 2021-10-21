import React, { FC } from 'react';
import { List, ListItem, ListItemText, Typography } from '@material-ui/core';

import { titleCase, isNumeric } from '../helpers';
import { LayerDefinition } from '../config/layers';

interface FeatureSidebarContentProps {
  f: any;
  layer: LayerDefinition;
}

export const FeatureSidebarContent: FC<FeatureSidebarContentProps> = ({ f, layer }) => {
  let details : React.ReactElement;
  switch (layer.id) {
    case 'airport_areas':
      details = <AirportDetails f={f} />;
      break;
    case 'elec_edges_high':
    case 'elec_edges_low':
      details = <PowerLineDetails f={f} />;
      break;
    case 'elec_nodes_junction':
      details = <PowerJunctionNodeDetails f={f} />;
      break;
    case 'elec_nodes_sink':
      details = <PowerDemandNodeDetails f={f} />;
      break;
    case 'elec_nodes_source':
      details = <PowerGenerationNodeDetails f={f} />;
      break;
    case 'port_areas':
      details = <PortDetails f={f} />;
      break;
    case 'rail_edges':
      details = <RailEdgeDetails f={f} />;
      break;
    case 'rail_nodes':
      details = <RailNodeDetails f={f} />;
      break;
    case 'road_bridges':
      details = <BridgeDetails f={f} />;
      break;
    case 'road_edges':
      details = <RoadEdgeDetails f={f} />;
      break;
    case 'water_irrigation_edges':
    case 'water_irrigation_nodes':
      details = <IrrigationDetails f={f} />;
      break;
    case 'water_potable_nodes':
      details = <WaterSupplyNodeDetails f={f} />;
      break;
    case 'water_potable_edges':
    case 'water_waste_edges':
      details = <WaterPipelineDetails f={f} />;
      break;
    case 'water_waste_nodes':
      details = <WastewaterNodeDetails f={f} />;
      break;
    default:
      details = <DefaultDetails f={f} />
  }

  return (
    <>
      <pre id="feature_debug" style={{ display: 'none' }}>
        <code>{JSON.stringify(f, null, 2)}</code>
      </pre>
      <Typography variant="caption">
        <span style={{ color: layer.color ?? '#333' }}>■</span>&nbsp;{layer.label}
      </Typography>
      {details}
    </>
  );
};

interface DataItemProps {
  label: string,
  value: any
}

const FORMATTER = Intl.NumberFormat('en-GB', { maximumSignificantDigits: 3 });

const DataItem: FC<DataItemProps> = ({ label, value }) => {
  if (isNumeric(value)) {
    value = FORMATTER.format(value)
  }
  return (
    <ListItem disableGutters>
      <ListItemText
        primary={label}
        primaryTypographyProps={{ variant: 'caption' }}
        secondary={value === ' ' ? '-' : value || '-'}
      />
    </ListItem>
  );
}

interface DetailSubheaderProps {
  id: string
}

const DetailSubheader: FC<DetailSubheaderProps> = ({ id }) => (
  <>
    <Typography variant="caption" component="p">
      ID: <span className="asset_id">{id}</span>
    </Typography>
  </>
);

interface DetailsProps {
  f: any
}

const DefaultDetails: FC<DetailsProps> = ({ f }) => {
  return (
    <>
      <Typography variant="h6" component="h1">Asset</Typography>
      <List>
        {Object.entries(f).map(([key, value]) => (
          <DataItem
            key={key}
            label={titleCase(key.replace(/_/g, ' '))}
            value={value}
          />
        ))}
      </List>
    </>
  );
}

const AirportDetails: FC<DetailsProps> = ({ f }) => (
  <>
    <Typography variant="h6" component="h1">{f.name}</Typography>
    <DetailSubheader id={f.asset_id} />
    <List>
      <DataItem
        label="Passengers (people/year)"
        value={f.passenger_number}
      />
      <DataItem
        label="Freight (tonnes/year)"
        value={f.freight_tonnes}
      />
      <DataItem
        label={`Rehabilitation cost (${f.unit})`}
        value={`${FORMATTER.format(f.maximum_damage)} (${FORMATTER.format(f.maximum_damage_lower)}–${FORMATTER.format(f.maximum_damage_upper)})`}
      />
      <DataItem
        label="Source"
        value={f.source}
      /></List>
  </>
);

const PowerLineDetails: FC<DetailsProps> = ({ f }) => (
  <>
    <Typography variant="h6" component="h1">{f.name ?? 'Power line'}</Typography>
    <DetailSubheader id={f.asset_id} />
    <List>
      <DataItem
        label="Connection"
        value={`${f.from_type} ${f.from_id}–${f.to_type} ${f.to_id} `}
      />
      <DataItem
        label="Voltage (kV)"
        value={f.voltage}
      />
      <DataItem
        label="Length (m)"
        value={f.length}
      />
      <DataItem
        label={`Rehabilitation cost (${f.cost_uom})`}
        value={`${FORMATTER.format(f.cost_avg)} (${FORMATTER.format(f.cost_min)}–${FORMATTER.format(f.cost_max)})`}
      />
      <DataItem
        label="Source"
        value={f.source}
      /></List>
  </>
);

const PowerGenerationNodeDetails: FC<DetailsProps> = ({ f }) => (
  <>
    <Typography variant="h6" component="h1">{`${f.title}–${f.subtype}`}</Typography>
    <DetailSubheader id={f.asset_id} />
    <List>
      <DataItem
        label="Capacity (MW)"
        value={f.capacity}
      />
      <DataItem
        label={`Energy intensity (${f.ei_uom})`}
        value={f.ei}
      />
      <DataItem
        label={`Rehabilitation cost (${f.cost_uom})`}
        value={`${FORMATTER.format(f.cost_avg)} (${FORMATTER.format(f.cost_min)}–${FORMATTER.format(f.cost_max)})`}
      />
      <DataItem
        label="Source"
        value={f.source}
      />
    </List>
  </>
);

const PowerDemandNodeDetails: FC<DetailsProps> = ({ f }) => (
  <>
    <Typography variant="h6" component="h1">{titleCase(f.subtype)}</Typography>
    <DetailSubheader id={f.asset_id} />
    <List>
      <DataItem
        label="Population served"
        value={f.population}
      />
      <DataItem
        label={`Energy intensity (${f.ei_uom})`}
        value={f.ei}
      />
      <DataItem
        label="Source"
        value={f.source}
      />
    </List>
  </>
);


const PowerJunctionNodeDetails: FC<DetailsProps> = ({ f }) => (
  <>
    <Typography variant="h6" component="h1">{titleCase(f.subtype)}</Typography>
    <DetailSubheader id={f.asset_id} />
    <List>
      <DataItem
        label={`Energy intensity (${f.ei_uom})`}
        value={f.ei}
      />
      <DataItem
        label={`Rehabilitation cost (${f.cost_uom})`}
        value={`${FORMATTER.format(f.cost_avg)} (${FORMATTER.format(f.cost_min)}–${FORMATTER.format(f.cost_max)})`}
      />
      <DataItem
        label="Source"
        value={f.source}
      />
    </List>
  </>
);

const IrrigationDetails: FC<DetailsProps> = ({ f }) => (
  <>
    <Typography variant="h6" component="h1">{f.asset}</Typography>
    <DetailSubheader id={f.asset_id} />
    <List>
      <DataItem
        label={`Rehabilitation cost (${f.cost_unit})`}
        value={`${FORMATTER.format(f.min_damage_cost)}–${FORMATTER.format(f.max_damage_cost)}`}
      />
      <DataItem
        label="Source"
        value={f.source}
      />
      <DataItem
        label="Notes"
        value={f.comment}
      />
    </List>
  </>
);

const WaterPipelineDetails: FC<DetailsProps> = ({ f }) => (
  <>
    <Typography variant="h6" component="h1">{f.asset}</Typography>
    <DetailSubheader id={f.asset_id} />
    <List>
      <DataItem
        label="Material"
        value={f.Material}
      />
      <DataItem
        label="Diameter (m)"
        value={f.Diameter}
      />
      <DataItem
        label="Length (m)"
        value={f.Length}
      />
      <DataItem
        label={`Rehabilitation cost (${f.cost_unit})`}
        value={`${FORMATTER.format(f.min_damage_cost)}–${FORMATTER.format(f.max_damage_cost)}`}
      />
      <DataItem
        label="Source"
        value={f.source}
      />
      <DataItem
        label="Notes"
        value={f.comment}
      />
    </List>
  </>
);

const PortDetails: FC<DetailsProps> = ({ f }) => (
  <>
    <Typography variant="h6" component="h1">{f.name}</Typography>
    <DetailSubheader id={f.asset_id} />
    <List>
      <DataItem
        label="Vessels (vessels/year)"
        value={f.vessels_number}
      />
      <DataItem
        label="Passengers (people/year)"
        value={f.passenger_number}
      />
      <DataItem
        label="Commodity"
        value={f.commodity}
      />
      <DataItem
        label="Export (tonnes/year)"
        value={f.export_tonnes}
      />
      <DataItem
        label="Import (tonnes/year)"
        value={f.import_tonnes}
      />
      <DataItem
        label="Import of vehicles (tonnes/year)"
        value={f.import_vehicles_tonnes}
      />
      <DataItem
        label="Transhipment (tonnes/year)"
        value={f.transhipment_tonnes}
      />
      <DataItem
        label={`Rehabilitation cost (${f.unit})`}
        value={`${FORMATTER.format(f.maximum_damage)} (${FORMATTER.format(f.maximum_damage_lower)}–${FORMATTER.format(f.maximum_damage_upper)})`}
      />
      <DataItem
        label="Source"
        value={f.source}
      />
      <DataItem
        label="Notes"
        value={f.comment}
      />
    </List>
  </>
);


const WaterSupplyNodeDetails: FC<DetailsProps> = ({ f }) => (
  <>
    <Typography variant="h6" component="h1">{f.name}</Typography>
    <DetailSubheader id={f.asset_id} />
    <List>
      <DataItem
        label="Vessels (vessels/year)"
        value={f.vessels_number}
      />
      <DataItem
        label="Tank"
        value={f.Tank_Type}
      />
      <DataItem
        label="Population served"
        value={f.asset_pop_new_2010}
      />
      <DataItem
        label="Capacity (mgd, millions of gallons/day)"
        value={f['capacity (mgd)']}
      />
      <DataItem
        label={`Rehabilitation cost (${f.cost_unit})`}
        value={`${FORMATTER.format(f.min_damage_cost)}–${FORMATTER.format(f.max_damage_cost)}`}
      />
      <DataItem
        label="Source"
        value={f.source}
      />
      <DataItem
        label="Notes"
        value={f.comment}
      />
    </List>
  </>
);

const RailEdgeDetails: FC<DetailsProps> = ({ f }) => (
  <>
    <Typography variant="h6" component="h1">{f.rail_sect}</Typography>
    <DetailSubheader id={f.asset_id} />
    <List>
      <DataItem
        label="Connection"
        value={`${f.from_node}–${f.to_node}`}
      />
      <DataItem
        label="Owner"
        value={f.type}
      />
      <DataItem
        label="Status"
        value={f.status}
      />
      <DataItem
        label="Length (m)"
        value={f.rail_length_m}
      />
      <DataItem
        label={`Rehabilitation cost (${f.cost_unit})`}value={`${FORMATTER.format(f.mean_damage_cost)} (${FORMATTER.format(f.min_damage_cost)}–${FORMATTER.format(f.max_damage_cost)})`}
      />
    </List>
  </>
);

const RailNodeDetails: FC<DetailsProps> = ({ f }) => (
  <>
    <Typography variant="h6" component="h1">{f.Station}</Typography>
    <DetailSubheader id={f.asset_id} />
    <List>
      <DataItem
        label="Lines"
        value={f.Lines}
      />
      <DataItem
        label="Status"
        value={`${f.status} (${f.Condition}) `}
      />
      <DataItem
        label={`Rehabilitation cost (${f.cost_unit})`}
        value={`${FORMATTER.format(f.mean_damage_cost)} (${FORMATTER.format(f.min_damage_cost)}–${FORMATTER.format(f.max_damage_cost)})`}
      />
    </List>
  </>
);

const RoadEdgeDetails: FC<DetailsProps> = ({ f }) => (
  <>
    <Typography variant="h6" component="h1">{f.street_name}{f.street_name? ', ' : ''}{f.section_name}</Typography>
    <DetailSubheader id={f.asset_id} />
    <List>
      <DataItem
        label="Connection"
        value={`${f.from_node}–${f.to_node}`}
      />
      <DataItem
        label="Road class (street type)"
        value={`${f.road_class} (${f.street_type ? f.street_type : 'none'})`}
      />
      <DataItem
        label="Construction"
        value={f.road_construction}
      />
      <DataItem
        label="Length (m)"
        value={f.road_length_m}
      />
      <DataItem
        label="Width (m)"
        value={f.road_width}
      />
      <DataItem
        label="Vertical alignment"
        value={f.vertalignm}
      />
      <DataItem
        label="Traffic (vehicles/day)"
        value={f.traffic_count}
      />
      <DataItem
        label={`Rehabilitation cost (${f.cost_unit})`}
        value={`${FORMATTER.format(f.mean_damage_cost)} (${FORMATTER.format(f.min_damage_cost)}–${FORMATTER.format(f.max_damage_cost)})`}
      />
      <DataItem
        label={`Reopening cost (${f.cost_unit})`}
        value={`${FORMATTER.format(f.mean_reopen_cost)} (${FORMATTER.format(f.min_reopen_cost)}–${FORMATTER.format(f.max_reopen_cost)})`}
      />
    </List>
  </>
);

const BridgeDetails: FC<DetailsProps> = ({ f }) => (
  <>
    <Typography variant="h6" component="h1">{titleCase(f.asset_type)}</Typography>
    <DetailSubheader id={f.asset_id} />
    <List>
      <DataItem
        label="Source ID"
        value={f.BRIDGEID}
      />
      <DataItem
        label={`Rehabilitation cost (${f.cost_unit})`}
        value={`${FORMATTER.format(f.mean_damage_cost)} (${FORMATTER.format(f.min_damage_cost)}–${FORMATTER.format(f.max_damage_cost)})`}
      />
      <DataItem
        label={`Reopening cost (${f.cost_unit})`}
        value={`${FORMATTER.format(f.mean_reopen_cost)} (${FORMATTER.format(f.min_reopen_cost)}–${FORMATTER.format(f.max_reopen_cost)})`}
      />
    </List>
  </>
);

const WastewaterNodeDetails: FC<DetailsProps> = ({ f }) => (
  <>
    <Typography variant="h6" component="h1">{f.Name}</Typography>
    <DetailSubheader id={f.asset_id} />
    <List>
      <DataItem
        label="Treatment type"
        value={f.Type_of_Tr}
      />
      <DataItem
        label="Status"
        value={f.Status}
      />
      <DataItem
        label="Capacity (mgd, millions of gallons/day)"
        value={f['capacity (mgd)']}
      />
      <DataItem
        label={`Rehabilitation cost (${f.cost_unit})`}
        value={`${FORMATTER.format(f.min_damage_cost)}–${FORMATTER.format(f.max_damage_cost)}`}
      />
      <DataItem
        label="Source"
        value={f.source}
      />
      <DataItem
        label="Notes"
        value={f.comment}
      />
    </List>
  </>
);
