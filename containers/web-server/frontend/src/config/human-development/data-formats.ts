import { FormatConfig } from '@/lib/data-map/view-layers';
import { toLabelLookup } from '@/lib/helpers';

import { HDI_VARIABLE_LABELS } from './metadata';

const hdiLabelLookup = toLabelLookup(HDI_VARIABLE_LABELS);

export function getHumanDevelopmentDataFormats(): FormatConfig {
  return {
    getDataLabel: ({ field }) => hdiLabelLookup[field],
    getValueFormatted: (value) => `${value}`,
  };
}
