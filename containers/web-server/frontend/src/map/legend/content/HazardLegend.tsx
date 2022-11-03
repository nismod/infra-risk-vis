import { FC, useCallback } from 'react';
import reactStringReplace from 'react-string-replace';

import { ViewLayer } from '@/lib/data-map/view-layers';

import { HAZARDS_METADATA, HAZARD_COLOR_MAPS, HazardType } from '@/config/hazards/metadata';

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

export const HazardLegend: FC<{ viewLayer: ViewLayer<{ hazardType: HazardType }> }> = ({ viewLayer }) => {
  const {
    params: { hazardType },
  } = viewLayer;

  let { label, dataUnit, labelAbbreviations = {}, legendAnnotation } = HAZARDS_METADATA[hazardType];
  const colorMap = HAZARD_COLOR_MAPS[hazardType];

  label = formatAbbreviations(label, labelAbbreviations);

  const getValueLabel = useCallback((value: number) => `${value.toLocaleString()} ${dataUnit}`, [dataUnit]);

  return (
    <RasterLegend label={label} description={legendAnnotation} colorMap={colorMap} getValueLabel={getValueLabel} />
  );
};
