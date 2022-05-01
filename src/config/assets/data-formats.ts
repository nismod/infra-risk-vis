import { HAZARDS_METADATA } from 'config/hazards/metadata';
import { FieldSpec, FormatConfig } from 'lib/data-map/view-layers';
import { numFormat, numFormatMoney, paren } from 'lib/helpers';
import _ from 'lodash';

function getSourceLabel(eadSource: string) {
  if (eadSource === 'all') return 'All Hazards';

  return HAZARDS_METADATA[eadSource].label;
}

const adaptationFieldLabels = {
  avoided_ead_mean: 'Avoided Direct Damages',
  avoided_eael_mean: 'Avoided Econ. Losses',
  adaptation_cost: 'Adaptation Cost',
  cost_benefit_ratio: 'Cost Benefit Ratio',
};

function getAdaptationFieldLabel(field: string) {
  return adaptationFieldLabels[field] || _.startCase(field);
}

const DATA_FORMATS: Record<string, FormatConfig> = {
  damages_expected: {
    getDataLabel: (colorField) => {
      const variableLabel = colorField.field === 'ead_mean' ? 'Direct Damages' : 'Economic Losses';
      const sourceLabel = getSourceLabel(colorField.fieldDimensions.hazard);
      return `${variableLabel} (${sourceLabel})`;
    },
    getValueFormatted: (value, fieldSpec) => {
      return value == null ? value : `${numFormatMoney(value)} $`;
    },
  },
  adaptation: {
    getDataLabel: (colorField) => {
      return `${getAdaptationFieldLabel(colorField.field)} ${paren(colorField.fieldDimensions.adaptation_name)}`;
    },
    getValueFormatted: (value, { field }) => {
      return value == null
        ? value
        : field === 'cost_benefit_ratio'
        ? `${numFormat(value)}x`
        : `${numFormatMoney(value)} $`;
    },
  },
};

export function getAssetDataFormats(fieldSpec: FieldSpec) {
  return DATA_FORMATS[fieldSpec.fieldGroup];
}
