import { Typography } from '@mui/material';
import { FC } from 'react';

import { InteractionTarget } from '@/lib/data-map/interactions/use-interactions';
import { numFormat } from '@/lib/helpers';
import { DataItem } from '@/lib/ui/data-display/DataItem';

import { REGIONS_METADATA } from '@/config/regions/metadata';

export const RegionDetailsContent: FC<{ selectedRegion: InteractionTarget<any> }> = ({ selectedRegion }) => {
  const metadata = REGIONS_METADATA[selectedRegion.viewLayer.params.regionLevel];
  const f = selectedRegion.target.feature.properties;
  const area = f['AREA'] ? f['AREA'] * 1e-6 : f['Shape_Area'] * 1e-6;

  return (
    <>
      {metadata.labelSingular}
      <Typography variant="h6">{selectedRegion.target.feature.properties[metadata.fieldName]}</Typography>
      <DataItem label="Population" value={f['population'].toLocaleString()} />
      <DataItem label="Area (km²)" value={numFormat(area)} />
      <DataItem label="Population per km²" value={numFormat(f['population_density_per_km2'])} />
    </>
  );
};
