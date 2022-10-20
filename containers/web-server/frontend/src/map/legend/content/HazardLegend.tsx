import { useCallback } from 'react';

import { RASTER_COLOR_MAPS } from '@/config/color-maps';
import { HAZARDS_METADATA } from '@/config/hazards/metadata';

import { RasterLegend } from '../RasterLegend';

export const HazardLegend = ({ viewLayer }) => {
  const {
    params: { hazardType },
  } = viewLayer;
  const { label, dataUnit } = HAZARDS_METADATA[hazardType];
  const colorMap = RASTER_COLOR_MAPS[hazardType];

  const getValueLabel = useCallback((value: number) => `${value.toLocaleString()} ${dataUnit}`, [dataUnit]);

  return <RasterLegend label={label} colorMap={colorMap} getValueLabel={getValueLabel} />;
};
