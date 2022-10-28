import { Typography } from '@mui/material';
import _ from 'lodash';
import { FC } from 'react';

import { InteractionTarget, VectorTarget } from '@/lib/data-map/interactions/use-interactions';
import { ColorBox } from '@/lib/ui/data-display/ColorBox';
import { DataDescription } from '@/lib/ui/data-display/DataDescription';
import { DataItem } from '@/lib/ui/data-display/DataItem';

import { TERRESTRIAL_COLORMAPS } from '@/config/_old/solutions/colors';
import { MARINE_HABITATS_LOOKUP } from '@/config/_old/solutions/domains';
import { habitatColorMap } from '@/state/layers/data-layers/_old/marine';
import { landuseColorMap } from '@/state/layers/data-layers/_old/terrestrial';

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
              colorSpec: TERRESTRIAL_COLORMAPS.slope,
            }}
            feature={feature}
            viewLayer={viewLayer}
          />
          <DataDescription
            colorMap={{
              fieldSpec: elevationFieldSpec,
              colorSpec: TERRESTRIAL_COLORMAPS.elevation,
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
