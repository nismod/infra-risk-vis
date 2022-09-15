import { FormatConfig } from 'lib/data-map/view-layers';
import { numFormat, numFormatMoney, toDictionary } from 'lib/helpers';
import {
  DroughtOptionsVariableType,
  DroughtRiskVariableType,
  DROUGHT_OPTIONS_VARIABLE_LABELS,
  DROUGHT_RISK_VARIABLE_LABELS,
} from './metadata';

const riskLabelLookup = toDictionary(
  DROUGHT_RISK_VARIABLE_LABELS,
  (x) => x.value,
  (x) => x.label,
);

const riskValueFormatLookup: Record<DroughtRiskVariableType, (x: any) => string> = {
  mean_monthly_water_stress_: (x) => numFormat(x),
  epd: (x) => numFormat(x),
  eael: (x) => `$${numFormat(x)}`,
};

export function getDroughtRiskDataFormats(): FormatConfig {
  return {
    getDataLabel: (colorField) => riskLabelLookup[colorField.field as DroughtRiskVariableType],
    getValueFormatted: (value, { field }) => riskValueFormatLookup[field](value),
  };
}

const optionsLabelLookup = toDictionary(
  DROUGHT_OPTIONS_VARIABLE_LABELS,
  (x) => x.value,
  (x) => x.label,
);

const optionsValueFormatLookup: Record<DroughtOptionsVariableType, (x: any) => string> = {
  cost_jmd: (x) => `$${numFormat(x)}`,
  population_protected: (x) => numFormat(x, 21),
  net_present_value_benefit: (x) => `$${numFormat(x)}`,
  benefit_cost_ratio: (x) => `${numFormat(x)}x`,
};

export function getDroughtOptionsDataFormats(): FormatConfig {
  return {
    getDataLabel: (colorField) => optionsLabelLookup[colorField.field as DroughtRiskVariableType],
    getValueFormatted: (value, { field }) => optionsValueFormatLookup[field](value),
  };
}
