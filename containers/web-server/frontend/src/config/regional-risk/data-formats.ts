import { FormatConfig } from '@/lib/data-map/view-layers';
import { toLabelLookup } from '@/lib/helpers';

import { REGIONAL_EXPOSURE_VARIABLE_LABELS } from './metadata';

const rexpLabelLookup = toLabelLookup(REGIONAL_EXPOSURE_VARIABLE_LABELS);

export function getRegionalExposureDataFormats(): FormatConfig {
  return {
    getDataLabel: ({ field }) => rexpLabelLookup[field],
    getValueFormatted: (value) => `${value}`,
  };
}
