import { Box, Typography } from '@mui/material';
import { VECTOR_COLOR_MAPS } from 'config/color-maps';
import { HAZARDS_METADATA } from 'config/hazards/metadata';
import { NETWORKS_METADATA } from 'config/networks/metadata';
import { DataItem } from 'features/detail-components';
import { colorMap } from 'lib/color-map';
import { InteractionTarget, VectorTarget } from 'lib/data-map/interactions/use-interactions';
import { ViewLayer } from 'lib/data-map/view-layers';
import { FC } from 'react';
import { useRecoilValue } from 'recoil';
import { damageSourceState, showDirectDamagesState } from 'state/damage-mapping/damage-map';
import { eadFieldSpecState } from 'state/damage-mapping/damage-style-params';
import { ColorBox } from './ColorBox';

function getSourceLabel(eadSource: string) {
  if (eadSource === 'total-damages') return 'Total';

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

const DirectDamagesDescription: FC<{ viewLayer: ViewLayer; feature: any }> = ({ viewLayer, feature }) => {
  const damageSource = useRecoilValue(damageSourceState);
  const eadFieldSpec = useRecoilValue(eadFieldSpecState);
  const eadAccessor = viewLayer.dataManager.getDataAccessor(viewLayer.id, eadFieldSpec);

  const value = eadAccessor?.(feature);
  const color = damageColorFn(value);
  return (
    <Box>
      <DataItem
        label={`Direct Damages (${getSourceLabel(damageSource)})`}
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
  const showDirectDamages = useRecoilValue(showDirectDamagesState);

  const f = hoveredObject.target.feature;
  const { label: title, color = '#ccc' } = NETWORKS_METADATA[hoveredObject.viewLayer.params.assetId];

  return (
    <>
      <Typography variant="body2">
        <ColorBox color={color} empty={showDirectDamages} />
        {title}
      </Typography>

      <DataItem label="ID" value={f.properties.asset_id} />
      {showDirectDamages && <DirectDamagesDescription viewLayer={hoveredObject.viewLayer} feature={f} />}
    </>
  );
};
