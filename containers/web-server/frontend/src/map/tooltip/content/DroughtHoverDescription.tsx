import { Typography } from '@mui/material';
import { FC, useMemo } from 'react';
import { useRecoilValue } from 'recoil';

import { InteractionTarget, VectorTarget } from '@/lib/data-map/interactions/use-interactions';

import { DataItem } from '@/details/features/detail-components';
import {
  droughtOptionsColorSpecState,
  droughtOptionsFieldSpecState,
  droughtRegionsColorSpecState,
  droughtRegionsFieldSpecState,
} from '@/state/layers/drought';

import { DataDescription } from '../DataDescription';

const DroughtRiskDescription: FC<{
  hoveredObject: InteractionTarget<VectorTarget>;
}> = ({ hoveredObject }) => {
  const {
    viewLayer,
    target: { feature },
  } = hoveredObject;

  const fieldSpec = useRecoilValue(droughtRegionsFieldSpecState);
  const colorSpec = useRecoilValue(droughtRegionsColorSpecState);

  const colorMap = useMemo(
    () => ({
      fieldSpec,
      colorSpec,
    }),
    [fieldSpec, colorSpec],
  );

  return (
    <>
      <Typography variant="body2">Drought Risk</Typography>

      <DataItem label="Region" value={feature.properties.HYDROLOGIC} />
      <DataDescription viewLayer={viewLayer} feature={feature} colorMap={colorMap} />
    </>
  );
};

const DroughtOptionDescription: FC<{
  hoveredObject: InteractionTarget<VectorTarget>;
}> = ({ hoveredObject }) => {
  const {
    viewLayer,
    target: { feature },
  } = hoveredObject;

  const fieldSpec = useRecoilValue(droughtOptionsFieldSpecState);
  const colorSpec = useRecoilValue(droughtOptionsColorSpecState);

  const colorMap = useMemo(
    () => ({
      fieldSpec,
      colorSpec,
    }),
    [fieldSpec, colorSpec],
  );

  return (
    <>
      <Typography variant="body2">Drought Adaptation Option</Typography>

      <DataItem label="Name" value={feature.properties.project_name} />
      <DataItem label="Type" value={feature.properties.project_type} />
      <DataDescription viewLayer={viewLayer} feature={feature} colorMap={colorMap} />
    </>
  );
};

export const DroughtHoverDescription: FC<{
  hoveredObject: InteractionTarget<VectorTarget>;
}> = ({ hoveredObject }) => {
  const { viewLayer } = hoveredObject;

  if (viewLayer.id === 'drought_risk') {
    return <DroughtRiskDescription hoveredObject={hoveredObject} />;
  } else if (viewLayer.id === 'drought_options') {
    return <DroughtOptionDescription hoveredObject={hoveredObject} />;
  }
};
