import { Box, Typography } from '@mui/material';
import { VECTOR_COLOR_MAPS } from 'config/color-maps';
import { HAZARDS_METADATA } from 'config/hazards/metadata';
import { NETWORKS_METADATA } from 'config/networks/metadata';
import { DataItem } from 'details/features/detail-components';
import { colorMap } from 'lib/color-map';
import { InteractionTarget, VectorTarget } from 'lib/data-map/interactions/use-interactions';
import { ViewLayer } from 'lib/data-map/view-layers';
import { FC, useMemo } from 'react';
import { useRecoilValue } from 'recoil';
import { damageSourceState, showDamagesState } from 'state/damage-mapping/damage-map';
import { damageMapStyleParamsState, damagesFieldState } from 'state/damage-mapping/damage-style-params';
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

const DamagesDescription: FC<{ viewLayer: ViewLayer; feature: any }> = ({ viewLayer, feature }) => {
  const damageSource = useRecoilValue(damageSourceState);
  const damagesFieldSpec = useRecoilValue(damagesFieldState);
  const styleParams = useRecoilValue(damageMapStyleParamsState);
  const eadAccessor = useMemo(() => viewLayer.dataAccessFn?.({ styleParams }).dataAccessor, [viewLayer, styleParams]);

  const value = eadAccessor?.(feature);
  const color = damageColorFn(value);
  const variableLabel = damagesFieldSpec.field === 'ead_mean' ? 'Direct Damages' : 'Economic Losses';
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

export const VectorHoverDescription: FC<{
  hoveredObject: InteractionTarget<VectorTarget>;
}> = ({ hoveredObject }) => {
  const showDirectDamages = useRecoilValue(showDamagesState);

  const {
    viewLayer,
    target: { feature },
  } = hoveredObject;

  const { label: title, color = '#ccc' } = NETWORKS_METADATA[viewLayer.params.assetId];

  return (
    <>
      <Typography variant="body2">
        <ColorBox color={color} empty={showDirectDamages} />
        {title}
      </Typography>

      <DataItem label="ID" value={feature.properties.asset_id} />
      {viewLayer.group === 'networks' && showDirectDamages && (
        <DamagesDescription viewLayer={viewLayer} feature={feature} />
      )}
    </>
  );
};
