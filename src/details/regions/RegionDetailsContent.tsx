import { Typography } from '@mui/material';
import { REGIONS_METADATA } from 'config/regions/metadata';
import { InteractionTarget } from 'lib/data-map/interactions/use-interactions';
import { FC } from 'react';

export const RegionDetailsContent: FC<{ selectedRegion: InteractionTarget<any> }> = ({ selectedRegion }) => {
  const metadata = REGIONS_METADATA[selectedRegion.viewLayer.params.boundaryLevel];

  return (
    <>
      <Typography variant="h6">{metadata.label}</Typography>
      {selectedRegion.target.feature.properties[metadata.fieldName]}
    </>
  );
};
