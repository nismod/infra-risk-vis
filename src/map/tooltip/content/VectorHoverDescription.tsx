import { Box, Typography } from '@mui/material';
import { HAZARDS_METADATA } from 'config/hazards/metadata';
import { NETWORKS_METADATA } from 'config/networks/metadata';
import { DataItem } from 'details/features/detail-components';
import { colorMap } from 'lib/color-map';
import { InteractionTarget, VectorTarget } from 'lib/data-map/interactions/use-interactions';
import { ColorMap, FieldSpec, ViewLayer } from 'lib/data-map/view-layers';
import { numFormat, paren } from 'lib/helpers';
import _ from 'lodash';
import { FC, useMemo } from 'react';
import { useRecoilValue } from 'recoil';
import { singleViewLayerParamsState } from 'state/layers/view-layers-params';
import { ColorBox } from './ColorBox';

function getSourceLabel(eadSource: string) {
  if (eadSource === 'all') return 'All Hazards';

  return HAZARDS_METADATA[eadSource].label;
}

function numFormatMoney(value: number) {
  return value.toLocaleString(undefined, {
    minimumFractionDigits: 2,
    maximumFractionDigits: 2,
  });
}

interface FormatConfig {
  getDataLabel: (viewLayer: ViewLayer, fieldSpec: FieldSpec) => string;
  getValueFormatted: (value: any, viewLayer: ViewLayer, fieldSpec: FieldSpec) => string;
}

const DATA_FORMATS: Record<string, FormatConfig> = {
  damages_expected: {
    getDataLabel: (viewLayer, colorField) => {
      const variableLabel = colorField.field === 'ead_mean' ? 'Direct Damages' : 'Economic Losses';
      const sourceLabel = getSourceLabel(colorField.fieldDimensions.hazard);
      return `${variableLabel} (${sourceLabel})`;
    },
    getValueFormatted: (value, viewLayer, fieldSpec) => {
      return value == null ? value : `${numFormatMoney(value)} $`;
    },
  },
  adaptation: {
    getDataLabel: (viewLayer, colorField) => {
      return `${_.startCase(colorField.field)} ${paren(colorField.fieldDimensions.adaptation_name)}`;
    },
    getValueFormatted: (value, viewLayer, { field }) => {
      return value == null
        ? value
        : field === 'cost_benefit_ratio'
        ? `${numFormat(value)}x`
        : `${numFormatMoney(value)} $`;
    },
  },
};

export const DataDescription: FC<{
  viewLayer: ViewLayer;
  feature: any;
  colorMap: ColorMap;
}> = ({ viewLayer, feature, colorMap: { fieldSpec: colorField, colorSpec } }) => {
  const accessor = useMemo(() => viewLayer.dataAccessFn?.(colorField), [viewLayer, colorField]);

  const value = accessor?.(feature);

  const colorFn = useMemo(() => colorMap(colorSpec), [colorSpec]);

  const color = colorFn(value);

  const { getDataLabel, getValueFormatted } = DATA_FORMATS[colorField.fieldGroup];

  const dataLabel = getDataLabel(viewLayer, colorField);
  const formattedValue = getValueFormatted(value, viewLayer, colorField);

  return (
    <Box>
      <DataItem
        label={dataLabel}
        value={
          <>
            <ColorBox color={color} />
            {formattedValue ?? '-'}
          </>
        }
      />
    </Box>
  );
};

export const VectorHoverDescription: FC<{
  hoveredObject: InteractionTarget<VectorTarget>;
}> = ({ hoveredObject }) => {
  const {
    viewLayer,
    target: { feature },
  } = hoveredObject;

  const layerParams = useRecoilValue(singleViewLayerParamsState(viewLayer.id));
  const { styleParams } = layerParams;
  const { colorMap } = styleParams ?? {};

  const isDataMapped = colorMap != null;

  const { label: title, color = '#ccc' } = NETWORKS_METADATA[viewLayer.params.assetId];

  return (
    <>
      <Typography variant="body2">
        <ColorBox color={color} empty={isDataMapped} />
        {title}
      </Typography>

      <DataItem label="ID" value={feature.properties.asset_id} />
      {colorMap && <DataDescription viewLayer={viewLayer} feature={feature} colorMap={colorMap} />}
    </>
  );
};
