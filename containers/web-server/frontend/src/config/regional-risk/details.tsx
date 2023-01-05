import { List, Typography } from '@mui/material';
import { FC } from 'react';

import { DataItem } from '@/lib/ui/data-display/DataItem';

import {
  DetailHeader,
  DetailsComponentProps,
  IdSubheader,
} from '@/details/features/detail-components';

import { REGIONAL_EXPOSURE_VARIABLE_LABELS } from './metadata';

export const RegionalExposureDetails: FC<DetailsComponentProps> = ({ f }) => {
  return (
    <>
      <DetailHeader>{f.NAME_EN}</DetailHeader>
      <IdSubheader id={f.ISO_A3} />
      <List>
        <Typography variant="subtitle2">Population exposed to:</Typography>
        {REGIONAL_EXPOSURE_VARIABLE_LABELS.map(({ value, label }) => (
          <DataItem label={label} value={f[value]} />
        ))}
      </List>
    </>
  );
};
