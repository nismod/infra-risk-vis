import { List } from '@mui/material';
import { FC } from 'react';

import { DataItem } from '@/lib/ui/data-display/DataItem';

import { DetailHeader, DetailsComponentProps, IdSubheader } from '@/details/features/detail-components';

export const HealthsiteDetails: FC<DetailsComponentProps> = ({ f }) => (
  <>
    <DetailHeader>{f.name}</DetailHeader>
    <IdSubheader id={f.osm_id} />
    <List>
      <DataItem label="Amenity" value={f.amenity} />
    </List>
  </>
);
