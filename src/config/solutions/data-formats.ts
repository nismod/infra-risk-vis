import { FieldSpec, FormatConfig } from 'lib/data-map/view-layers';

const variableLabelLookup = {
  landuse_desc: 'Land Use',
  slope_degrees: 'Slope (deg)',
  elevation_m: 'Elevation (m)',
};

const valueFormatLookup = {
  slope_degrees: (x) => x.toLocaleString(undefined, { maximumFractionDigits: 0 }),
  elevation_m: (x) => x,
};

export function getTerrestrialDataFormats(fieldSpec: FieldSpec): FormatConfig {
  return {
    getDataLabel: ({ field }) => variableLabelLookup[field],
    getValueFormatted: (value, { field }) => valueFormatLookup[field](value),
  };
}
