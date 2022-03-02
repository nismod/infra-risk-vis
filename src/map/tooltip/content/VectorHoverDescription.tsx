import { Box, Typography } from '@mui/material';
import { VECTOR_COLOR_MAPS } from 'config/color-maps';
import { HAZARDS_METADATA } from 'config/hazards/metadata';
import { NETWORKS_METADATA } from 'config/networks/metadata';
import { DataItem } from 'features/detail-components';
import { InteractionTarget, VectorTarget } from 'lib/data-map/interactions/use-interactions';
import { colorMap } from 'lib/deck-layers/utils';
import { FC } from 'react';
import { useRecoilValue } from 'recoil';
import { selectedDamageSourceState, showDirectDamagesState } from 'state/damage-mapping/damage-map';
import { eadAccessorState } from 'state/damage-mapping/damage-style-params';
import { ColorBox } from './ColorBox';

function getSourceLabel(eadSource: string) {
  if (eadSource === 'total-damages') return 'Total';

  return HAZARDS_METADATA[eadSource].label;
}

function formatDamagesValue(value: number) {
  if (value == null) return '-';

  return value.toLocaleString(undefined, { minimumFractionDigits: 2, maximumFractionDigits: 2 }) + ' $';
}

const damageColorSpec = VECTOR_COLOR_MAPS['damages'];
const damageColorFn = colorMap(damageColorSpec.scale, damageColorSpec.range, damageColorSpec.empty);

const DirectDamagesDescription: FC<{ feature: any }> = ({ feature }) => {
  const damageSource = useRecoilValue(selectedDamageSourceState);
  const eadAccessor = useRecoilValue(eadAccessorState);

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

export const VectorHoverDescription: FC<{ hoveredObject: InteractionTarget<VectorTarget> }> = ({ hoveredObject }) => {
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
      {showDirectDamages && <DirectDamagesDescription feature={f} />}
    </>
  );
};
