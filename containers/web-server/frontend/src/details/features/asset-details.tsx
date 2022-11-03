import { Download } from '@mui/icons-material';
import { IconButton, Typography } from '@mui/material';
import { Box } from '@mui/system';
import { FC, ReactElement, Suspense } from 'react';
import { RecoilValue, useRecoilValue } from 'recoil';

import { downloadFile } from '@/lib/helpers';
import { ColorBox } from '@/lib/ui/data-display/ColorBox';

import { apiFeatureQuery } from '@/state/queries';

import { AdaptationSection } from './adaptation/AdaptationSection';
import { DamagesSection } from './damages/DamagesSection';
import { DetailsComponentType } from './detail-components';

interface LoadDetailsProps {
  featureDetailsState: RecoilValue<any>;
  children: (details: any) => ReactElement;
}

const LoadDetails: FC<LoadDetailsProps> = ({ featureDetailsState, children }) => {
  const featureDetails = useRecoilValue(featureDetailsState);

  return children(featureDetails);
};

const DownloadButton = ({ feature }) => {
  return (
    <IconButton
      sx={{
        position: 'absolute',
        top: 0,
        right: 30, // hack: larger right margin to allow space for close button
      }}
      title="Download CSV with feature metadata"
      onClick={() => downloadFile(makeDetailsCsv(feature), 'text/csv', `feature_${feature.id}.csv`)}
    >
      <Download />
    </IconButton>
  );
};

export interface AssetDetailsProps {
  feature: any;
  label: string;
  color?: string;
}

const AssetDetails: FC<AssetDetailsProps> = ({ label, color = '#cccccc', feature, children }) => {
  return (
    <Box position="relative">
      <code style={{ display: 'none' }} className="feature-debug">
        {JSON.stringify(feature.properties, null, 2)}
      </code>
      <Typography variant="caption">
        <ColorBox color={color} />
        {label}
      </Typography>
      <Box>{children}</Box>
    </Box>
  );
};

type SimpleAssetDetailsProps = AssetDetailsProps & {
  detailsComponent: DetailsComponentType;
};

export const SimpleAssetDetails: FC<SimpleAssetDetailsProps> = ({ label, color, detailsComponent, feature }) => {
  const DetailsComponent = detailsComponent;

  return (
    <AssetDetails label={label} color={color} feature={feature}>
      <DetailsComponent f={feature.properties} />
      <DownloadButton feature={feature} />
    </AssetDetails>
  );
};

type ExtendedAssetDetailsProps = SimpleAssetDetailsProps & {
  showRiskSection: boolean;
};

export const ExtendedAssetDetails: FC<ExtendedAssetDetailsProps> = ({
  label,
  color,
  detailsComponent,
  feature,
  showRiskSection,
}) => {
  const DetailsComponent = detailsComponent;

  const featureDetailsState = apiFeatureQuery(feature.id);

  return (
    <AssetDetails label={label} color={color} feature={feature}>
      <Suspense fallback="Loading data...">
        <LoadDetails featureDetailsState={featureDetailsState}>
          {(featureDetails) => (
            <>
              <code style={{ display: 'none' }} className="feature-details-debug">
                {JSON.stringify(featureDetails.properties, null, 2)}
              </code>
              <DetailsComponent f={featureDetails.properties} />
              <DownloadButton feature={featureDetails} />
              {showRiskSection && (
                <>
                  <DamagesSection fd={featureDetails} />
                  <AdaptationSection fd={featureDetails} />
                </>
              )}
            </>
          )}
        </LoadDetails>
      </Suspense>
    </AssetDetails>
  );
};

function makeDetailsCsv(fd) {
  return (
    'variable,value\n' +
    Object.entries(fd.properties)
      .map(([k, v]) => `${k},${v}`)
      .join('\n')
  );
}
