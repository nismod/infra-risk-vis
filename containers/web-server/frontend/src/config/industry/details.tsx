import { List } from '@mui/material';

import { DataItem } from '@/lib/ui/data-display/DataItem';

import { DetailHeader, IdSubheader } from '@/details/features/detail-components';

export const IndustryDetails = ({ f }) => (
  <>
    <DetailHeader>{f.name}</DetailHeader>
    <IdSubheader id={f.uid} />
    <List>
      <DataItem label="Country" value={f.country} />
      <DataItem label="Primary Product" value={f.primary_product} />
    </List>
  </>
);
