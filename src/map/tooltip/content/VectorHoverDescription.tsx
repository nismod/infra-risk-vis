import { Box, Typography } from '@mui/material';
import { VECTOR_COLOR_MAPS } from 'config/color-maps';
import { HAZARDS_METADATA } from 'config/hazards/metadata';
import { NETWORKS_METADATA } from 'config/networks/metadata';
import { DataItem } from 'details/features/detail-components';
import { colorMap } from 'lib/color-map';
import { InteractionTarget, VectorTarget } from 'lib/data-map/interactions/use-interactions';
import { StyleParams, ViewLayer } from 'lib/data-map/view-layers';
import { paren } from 'lib/helpers';
import _ from 'lodash';
import { FC, useMemo } from 'react';
import { useRecoilValue } from 'recoil';
import { singleViewLayerParamsState } from 'state/layers/view-layers-params';
import { ColorBox } from './ColorBox';

function getSourceLabel(eadSource: string) {
  if (eadSource === 'all') return 'All Hazards';

  return HAZARDS_METADATA[eadSource].label;
}

function formatDamagesValue(value: number) {
  if (value == null) return '-';

  return (
    value.toLocaleString(undefined, {
      minimumFractionDigits: 2,
      maximumFractionDigits: 2,
    }) + ' $'
  );
}

const damageColorSpec = VECTOR_COLOR_MAPS['damages'];
const damageColorFn = colorMap(damageColorSpec.scale, damageColorSpec.range, damageColorSpec.empty);

const DamagesDescription: FC<{ viewLayer: ViewLayer; feature: any; styleParams: StyleParams }> = ({
  viewLayer,
  feature,
  styleParams,
}) => {
  const eadAccessor = useMemo(() => viewLayer.dataAccessFn?.({ styleParams }).dataAccessor, [viewLayer, styleParams]);

  const fieldSpec = styleParams.colorMap.colorField;
  const damageSource = styleParams.colorMap.colorField.fieldDimensions.hazard;
  const value = eadAccessor?.(feature);
  const color = damageColorFn(value);
  const variableLabel = fieldSpec.field === 'ead_mean' ? 'Direct Damages' : 'Economic Losses';
  return (
    <Box>
      <DataItem
        label={`${variableLabel} (${getSourceLabel(damageSource)})`}
        value={
          <>
            <ColorBox color={color} />
            {formatDamagesValue(value)}
          </>
        }
      />
    </Box>
  );
};

const AdaptationDescription: FC<{ viewLayer: ViewLayer; feature: any; styleParams: StyleParams }> = ({
  viewLayer,
  feature,
  styleParams,
}) => {
  const {
    colorMap: { colorField, colorScheme },
  } = styleParams;
  const adaptationAccessor = useMemo(
    () => viewLayer.dataAccessFn?.({ styleParams }).dataAccessor,
    [viewLayer, styleParams],
  );
  const value = adaptationAccessor?.(feature);

  const colorSpec = VECTOR_COLOR_MAPS[colorScheme];
  const colorFn = useMemo(() => colorMap(colorSpec.scale, colorSpec.range, colorSpec.empty), [colorSpec]);

  const color = colorFn(value);

  return (
    <Box>
      <DataItem
        label={`${_.startCase(colorField.field)} ${paren(colorField.fieldDimensions.adaptation_name)}`}
        value={
          <>
            <ColorBox color={color} />
            {formatDamagesValue(value)}
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
  const { colorMap } = styleParams;

  const isDataMapped = colorMap != null;

  const { label: title, color = '#ccc' } = NETWORKS_METADATA[viewLayer.params.assetId];

  return (
    <>
      <Typography variant="body2">
        <ColorBox color={color} empty={isDataMapped} />
        {title}
      </Typography>

      <DataItem label="ID" value={feature.properties.asset_id} />
      {colorMap?.colorField.fieldGroup === 'damages_expected' && (
        <DamagesDescription viewLayer={viewLayer} feature={feature} styleParams={styleParams} />
      )}
      {colorMap?.colorField.fieldGroup === 'adaptation' && (
        <AdaptationDescription viewLayer={viewLayer} feature={feature} styleParams={styleParams} />
      )}
    </>
  );
};
