import { Typography } from '@mui/material';
import { VECTOR_COLOR_MAPS } from 'config/color-maps';
import { MARINE_HABITATS_LOOKUP } from 'config/solutions/domains';
import { DataItem } from 'details/features/detail-components';
import { colorMap } from 'lib/color-map';
import { InteractionTarget, VectorTarget } from 'lib/data-map/interactions/use-interactions';
import _ from 'lodash';
import { FC } from 'react';
import { habitatColorMap } from 'state/layers/marine';
import { landuseColorMap } from 'state/layers/terrestrial';
import { ColorBox } from './ColorBox';

const slopeColorFunction = colorMap(VECTOR_COLOR_MAPS.terrestrialSlope);
const elevationColorFunction = colorMap(VECTOR_COLOR_MAPS.terrestrialElevation);

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
          <DataItem
            label="Land Use"
            value={
              <>
                <ColorBox color={landuseColorMap(feature.properties.landuse_desc)} />
                {feature.properties.landuse_desc}
              </>
            }
          />
          <DataItem
            label="Slope (deg)"
            value={
              <>
                <ColorBox color={slopeColorFunction(feature.properties.slope_degrees)} />
                {feature.properties.slope_degrees.toLocaleString(undefined, { maximumFractionDigits: 0 })}
              </>
            }
          />
          <DataItem
            label="Elevation (m)"
            value={
              <>
                <ColorBox color={elevationColorFunction(feature.properties.elevation_m)} />
                {feature.properties.elevation_m}
              </>
            }
          />
        </>
      )}
      {viewLayer.id === 'marine' && (
        <>
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
