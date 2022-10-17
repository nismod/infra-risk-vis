import { Typography } from '@mui/material';
import _ from 'lodash';
import { FC } from 'react';

import { colorMap } from '@/lib/color-map';

import { VECTOR_COLOR_MAPS } from '@/config/color-maps';
import { MARINE_HABITATS_LOOKUP } from '@/config/solutions/domains';
import { DataItem } from '@/details/features/detail-components';
import { ColorBox } from '@/map/tooltip/content/ColorBox';
import { habitatColorMap } from '@/state/layers/marine';
import { landuseColorMap } from '@/state/layers/terrestrial';

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
          {feature.properties.is_coral ? <DataItem label="Coral Type" value={feature.properties.coral_type} /> : null}
          {feature.properties.is_mangrove ? (
            <DataItem label="Mangrove Type" value={feature.properties.mangrove_type} />
          ) : null}
          <Typography variant="caption" component="h6">
            Proximity
          </Typography>
          {feature.properties.within_coral_500m ? (
            <Typography variant="body2" component="p">
              within 500m of coral
            </Typography>
          ) : null}
          {feature.properties.within_mangrove_500m ? (
            <Typography variant="body2" component="p">
              within 500m of mangrove
            </Typography>
          ) : null}
          {feature.properties.within_seagrass_500m ? (
            <Typography variant="body2" component="p">
              within 500m of seagrass
            </Typography>
          ) : null}
        </>
      )}
    </>
  );
};
