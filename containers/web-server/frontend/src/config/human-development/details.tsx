import { List } from '@mui/material';
import { FC } from 'react';

import { DataItem } from '@/lib/ui/data-display/DataItem';

import {
  DetailHeader,
  DetailSubheader,
  DetailsComponentProps,
  DetailsComponentType,
} from '@/details/features/detail-components';

import { HdiRegionLevel } from './metadata';

const HDIDataList = ({ f }) => (
  <List>
    <DataItem label="Human Development Index" value={f.subnational_hdi} />
    <DataItem label="Educational Index" value={f.educational_index} />
    <DataItem label="Income Index" value={f.income_index} />
    <DataItem label="Health Index" value={f.health_index} />
  </List>
);

export const HdiCountryDetails: FC<DetailsComponentProps> = ({ f }) => {
  return (
    <>
      <DetailHeader>{f.country}</DetailHeader>
      <HDIDataList f={f} />
    </>
  );
};

export const HdiRegionDetails: FC<DetailsComponentProps> = ({ f }) => {
  return (
    <>
      <DetailHeader>{f.region}</DetailHeader>
      <DetailSubheader>Country: {f.country}</DetailSubheader>
      <HDIDataList f={f} />
    </>
  );
};

export const HDI_REGION_LEVEL_DETAILS: Record<HdiRegionLevel, DetailsComponentType> = {
  countries: HdiCountryDetails,
  regions: HdiRegionDetails,
};
