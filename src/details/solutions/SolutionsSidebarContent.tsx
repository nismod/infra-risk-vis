import { Typography } from '@mui/material';
import { VECTOR_COLOR_MAPS } from 'config/color-maps';
import { MARINE_HABITATS_LOOKUP } from 'config/solutions/domains';
import { DataItem } from 'details/features/detail-components';
import { colorMap } from 'lib/color-map';
import _ from 'lodash';
import { ColorBox } from 'map/tooltip/content/ColorBox';
import { FC } from 'react';
import { habitatColorMap } from 'state/layers/marine';
import { landuseColorMap } from 'state/layers/terrestrial';

const slopeColorFunction = colorMap(VECTOR_COLOR_MAPS.terrestrialSlope);
const elevationColorFunction = colorMap(VECTOR_COLOR_MAPS.terrestrialElevation);

interface SolutionsSidebarContentProps {
  feature: any;
  solutionType: string; // 'terrestrial' | 'marine'
  showRiskSection?: boolean;
}

export const SolutionsSidebarContent: FC<SolutionsSidebarContentProps> = ({ feature, solutionType }) => {
  return (
    <>
      <Typography variant="body2">{_.startCase(solutionType)}</Typography>

      {solutionType === 'terrestrial' && (
        <>
          <DataItem label="Cell ID" value={feature.properties.cell_index} maximumSignificantDigits={21} />
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
      {solutionType === 'marine' && (
        <DataItem
          label="Habitat"
          value={
            <>
              <ColorBox color={habitatColorMap(feature.properties.habitat)} />
              {feature.properties.habitat ? MARINE_HABITATS_LOOKUP[feature.properties.habitat] : 'Buffer Zone'}
            </>
          }
        />
      )}
    </>
  );
};
