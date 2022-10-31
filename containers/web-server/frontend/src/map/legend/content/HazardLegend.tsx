import { useCallback } from 'react';
import reactStringReplace from 'react-string-replace';

import { RASTER_COLOR_MAPS } from '@/config/color-maps';
import { HAZARDS_METADATA } from '@/config/hazards/metadata';

import { RasterLegend } from '../RasterLegend';

function formatAbbreviations(label, abbreviations) {
  for (let [key, value] of Object.entries(abbreviations)) {
    label = reactStringReplace(label, key, (match, i) => (
      <abbr key={match + i} title={value as string}>
        {key}
      </abbr>
    ));
  }

  return label;
}

export const HazardLegend = ({ viewLayer }) => {
  const {
    params: { hazardType },
  } = viewLayer;
  let { label, dataUnit, labelAbbreviations = {} } = HAZARDS_METADATA[hazardType];
  const colorMap = RASTER_COLOR_MAPS[hazardType];

  label = formatAbbreviations(label, labelAbbreviations);

  const getValueLabel = useCallback((value: number) => `${value.toLocaleString()} ${dataUnit}`, [dataUnit]);

  return <RasterLegend label={label} colorMap={colorMap} getValueLabel={getValueLabel} />;
};
