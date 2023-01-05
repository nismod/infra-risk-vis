import { Typography } from '@mui/material';
import { Box } from '@mui/system';
import { FC, ReactElement, Suspense } from 'react';
import { RecoilValue, useRecoilValue } from 'recoil';

import { ColorBox } from '@/lib/ui/data-display/ColorBox';

import { apiFeatureQuery } from '@/state/queries';

import { ButtonPlacement, DownloadButton } from './DownloadButton';
// import { AdaptationSection } from './adaptation/AdaptationSection';
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

export const SimpleAssetDetails: FC<SimpleAssetDetailsProps> = ({
  label,
  color,
  detailsComponent,
  feature,
}) => {
  const DetailsComponent = detailsComponent;

  return (
    <AssetDetails label={label} color={color} feature={feature}>
      <DetailsComponent f={feature.properties} />
      <ButtonPlacement
        right={30} //hack: larger right margin to allow space for close button
      >
        <DownloadButton
          makeContent={() => makeDetailsCsv(feature)}
          title="Download CSV with feature metadata"
          filename={`feature_${feature.id}.csv`}
        />
      </ButtonPlacement>
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
      <ButtonPlacement
        right={30} // hack: larger right margin to allow space for close button
      >
        <DownloadButton
          makeContent={() => makeDetailsCsv(feature)}
          title="Download CSV with feature metadata"
          filename={`feature_${feature.id}.csv`}
        />
      </ButtonPlacement>
      <Suspense fallback="Loading data...">
        <LoadDetails featureDetailsState={featureDetailsState}>
          {(featureDetails) => (
            <>
              <DetailsComponent f={featureDetails.properties} />
              {showRiskSection && (
                <>
                  <DamagesSection fd={featureDetails} />
                  {/* <AdaptationSection fd={featureDetails} /> */}
                </>
              )}
              <details className="feature-details-debug">
                <summary>
                  <small>Feature data</small>
                </summary>
                <pre>{JSON.stringify(featureDetails, null, 2)}</pre>
              </details>
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
