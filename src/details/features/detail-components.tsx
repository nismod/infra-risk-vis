import { ComponentType, FC } from 'react';

import { List, ListItem, ListItemText, Typography } from '@mui/material';

import { titleCase, isNumeric, numFormat, paren, numRangeFormat } from 'lib/helpers';

interface DataItemProps {
  label: string;
  value: any;
  maximumSignificantDigits?: number;
}

export const DataItem: FC<DataItemProps> = ({ label, value, maximumSignificantDigits }) => {
  if (isNumeric(value)) {
    value = numFormat(value, maximumSignificantDigits);
  }
  return (
    <ListItem disableGutters disablePadding>
      <ListItemText
        primary={label}
        primaryTypographyProps={{ variant: 'caption' }}
        secondary={value === ' ' ? '-' : value || '-'}
      />
    </ListItem>
  );
};

interface DetailSubheaderProps {
  id: string;
}

const DetailSubheader: FC<DetailSubheaderProps> = ({ id }) => (
  <Typography variant="caption" component="p">
    ID: <span className="asset_id">{id}</span>
  </Typography>
);

interface DetailsComponentProps {
  f: any;
}

export type DetailsComponent = ComponentType<DetailsComponentProps>;

export const DefaultDetailsList: FC<DetailsComponentProps> = ({ f }) => {
  return (
    <List>
      {Object.entries(f).map(([key, value]) => (
        <DataItem key={key} label={titleCase(key.replace(/_/g, ' '))} value={value} />
      ))}
    </List>
  );
};

export const DefaultDetails: FC<DetailsComponentProps> = ({ f }) => {
  return (
    <>
      <Typography variant="h6" component="h1">
        Asset
      </Typography>
      <DefaultDetailsList f={f} />
    </>
  );
};

export const AirportDetails: FC<DetailsComponentProps> = ({ f }) => (
  <>
    <Typography variant="h6" component="h1">
      {f.name}
    </Typography>
    <DetailSubheader id={f.asset_id} />
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
    {f.name && (
      <Typography variant="h6" component="h1">
        {f.name}
      </Typography>
    )}
    <DetailSubheader id={f.asset_id} />
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
    <Typography variant="h6" component="h1">
      {f.title}
    </Typography>
    <DetailSubheader id={f.asset_id} />
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
    <DetailSubheader id={f.asset_id} />
    <List>
      <DataItem label="Population served" value={f.population} />
      <DataItem label={`Energy intensity (${f.ei_uom})`} value={f.ei} />
    </List>
  </>
);

export const PowerJunctionNodeDetails: FC<DetailsComponentProps> = ({ f }) => (
  <>
    <DetailSubheader id={f.asset_id} />
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
    <DetailSubheader id={f.asset_id} />
    <List>
      <DataItem label={`Rehabilitation cost (${f.cost_unit})`} value={numRangeFormat(f.cost_min, f.cost_max)} />
      <DataItem label="Notes" value={f.comment} />
    </List>
  </>
);

export const WaterPipelineDetails: FC<DetailsComponentProps> = ({ f }) => (
  <>
    <DetailSubheader id={f.asset_id} />
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
    <Typography variant="h6" component="h1">
      {f.name}
    </Typography>
    <DetailSubheader id={f.asset_id} />
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
    <Typography variant="h6" component="h1">
      {f.name}
    </Typography>
    <DetailSubheader id={f.asset_id} />
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
    <Typography variant="h6" component="h1">
      {f.rail_sect}
    </Typography>
    <DetailSubheader id={f.asset_id} />
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
    <Typography variant="h6" component="h1">
      {f.Station}
    </Typography>
    <DetailSubheader id={f.asset_id} />
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
    <Typography variant="h6" component="h1">
      {f.street_name}
      {f.street_name ? ', ' : ''}
      {f.section_name}
    </Typography>
    <DetailSubheader id={f.asset_id} />
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
    <Typography variant="h6" component="h1">
      {f.name}
    </Typography>
    <DetailSubheader id={f.asset_id} />
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
    <DetailSubheader id={f.asset_id} />
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
    <Typography variant="h6" component="h1">
      {f.Name}
    </Typography>
    <DetailSubheader id={f.asset_id} />
    <List>
      <DataItem label="Treatment type" value={f.Type_of_Tr} />
      <DataItem label="Status" value={f.Status} />
      <DataItem label="Capacity (mgd, millions of gallons/day)" value={f['capacity (mgd)']} />
      <DataItem label={`Rehabilitation cost (${f.cost_unit})`} value={numRangeFormat(f.cost_min, f.cost_max)} />
      <DataItem label="Notes" value={f.comment} />
    </List>
  </>
);

export const BuildingDetails: FC<DetailsComponentProps> = ({ f }) => (
  <>
    {f.name && (
      <Typography variant="h6" component="h1">
        {f.name}
      </Typography>
    )}
    <DetailSubheader id={f.asset_id} />
    <List>
      <DataItem label="Source ID" value={f.osm_way_id} />
      <DataItem label={`Total GDP (${f.GDP_unit})`} value={numFormat(f.total_GDP)} />
      <DataItem
        label={`Rehabilitation cost (${f.cost_unit})`}
        value={`${numFormat(f.cost_mean)} ${paren(numRangeFormat(f.cost_min, f.cost_max))}`}
      />
    </List>
  </>
);
