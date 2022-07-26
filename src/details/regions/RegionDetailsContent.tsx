import { Typography } from '@mui/material';
import { REGIONS_METADATA } from 'config/regions/metadata';
import { DataItem } from 'details/features/detail-components';
import { InteractionTarget } from 'lib/data-map/interactions/use-interactions';
import { numFormat } from 'lib/helpers';
import { FC } from 'react';

export const RegionDetailsContent: FC<{ selectedRegion: InteractionTarget<any> }> = ({ selectedRegion }) => {
  const metadata = REGIONS_METADATA[selectedRegion.viewLayer.params.regionLevel];

  return (
    <>
      {metadata.labelSingular}
      <Typography variant="h6">{selectedRegion.target.feature.properties[metadata.fieldName]}</Typography>
    </>
  );
};
