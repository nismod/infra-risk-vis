import { Typography } from '@mui/material';
import _ from 'lodash';
import { FC } from 'react';

import { InteractionTarget, VectorTarget } from '@/lib/data-map/interactions/use-interactions';

import { VECTOR_COLOR_MAPS } from '@/config/color-maps';
import { MARINE_HABITATS_LOOKUP } from '@/config/solutions/domains';
import { DataItem } from '@/details/features/detail-components';
import { habitatColorMap } from '@/state/layers/marine';
import { landuseColorMap } from '@/state/layers/terrestrial';

import { DataDescription } from '../DataDescription';
import { ColorBox } from './ColorBox';

const slopeFieldSpec = {
  fieldGroup: 'properties',
  field: 'slope_degrees',
};

const elevationFieldSpec = {
  fieldGroup: 'properties',
  field: 'elevation_m',
};

export const SolutionHoverDescription: FC<{
  hoveredObject: InteractionTarget<VectorTarget>;
}> = ({ hoveredObject }) => {
  const {
    viewLayer,
    target: { feature },
  } = hoveredObject;

  return (
    <>
      <Typography variant="body2">{_.startCase(viewLayer.id)}</Typography>

      {viewLayer.id === 'terrestrial' && (
        <>
          <DataItem label="Cell ID" value={feature.properties.cell_id} maximumSignificantDigits={21} />
          {/* not using DataDescription for Land Use because currently it only works for colorSpec-based color maps (not categorical) */}
          <DataItem
            label="Land Use"
            value={
              <>
                <ColorBox color={landuseColorMap(feature.properties.landuse_desc)} />
                {feature.properties.landuse_desc}
              </>
            }
          />
          <DataDescription
            colorMap={{
              fieldSpec: slopeFieldSpec,
              colorSpec: VECTOR_COLOR_MAPS.terrestrialSlope,
            }}
            feature={feature}
            viewLayer={viewLayer}
          />
          <DataDescription
            colorMap={{
              fieldSpec: elevationFieldSpec,
              colorSpec: VECTOR_COLOR_MAPS.terrestrialElevation,
            }}
            feature={feature}
            viewLayer={viewLayer}
          />
        </>
      )}
      {viewLayer.id === 'marine' && (
        <>
          {/* not using DataDescription for Habitat because currently it only works for colorSpec-based color maps (not categorical) */}
          <DataItem
            label="Habitat"
            value={
              <>
                <ColorBox color={habitatColorMap(feature.properties.habitat)} />
                {feature.properties.habitat ? MARINE_HABITATS_LOOKUP[feature.properties.habitat] : 'Buffer Zone'}
              </>
            }
          />
        </>
      )}
    </>
  );
};
