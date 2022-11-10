import { Box, Fade, Tooltip, Typography } from '@mui/material';
import { FC, ReactNode } from 'react';

import { RasterColorMapValues } from './RasterLegend';

const legendHeight = 10;

export function formatRangeTruncation(
  value: ReactNode | string,
  i: number,
  [minTruncated, maxTruncated]: [boolean, boolean],
) {
  if (i === 0 && minTruncated) {
    return <>&le;{value}</>;
  }
  // terracotta sends 255 color values, so last index is 254
  if (i === 254 && maxTruncated) {
    return <>&ge;{value}</>;
  }
  return value;
}

const LegendGradient: FC<{
  colorMap: RasterColorMapValues;
  getValueLabel: (value: number) => ReactNode | string;
}> = ({ colorMap, getValueLabel }) => {
  const { colorMapValues, rangeTruncated } = colorMap;
  return (
    <>
      {colorMapValues.map(({ color, value }, i) => (
        <Tooltip
          key={i}
          title={formatRangeTruncation(getValueLabel(value), i, rangeTruncated)}
          arrow
          placement="top"
          enterDelay={200}
          leaveDelay={0}
          // deactivates transition animation
          TransitionComponent={Fade}
          TransitionProps={{ timeout: 0 }}
        >
          <Box height={legendHeight} flexGrow={1} bgcolor={color} />
        </Tooltip>
      ))}
    </>
  );
};

export interface GradientLegendProps {
  label: string | ReactNode;
  description?: string;
  range: [number, number];
  colorMap: RasterColorMapValues;
  getValueLabel: (x: any) => ReactNode | string;
}

export const GradientLegend: FC<GradientLegendProps> = ({ label, description, range, colorMap, getValueLabel }) => (
  <Box mb={2}>
    <Box mb={1}>
      <Typography variant="body1">{label}</Typography>
      {description && <Typography variant="body2">{description}</Typography>}
    </Box>
    <Box
      height={legendHeight + 2}
      width={256}
      bgcolor="#ccc"
      display="flex"
      flexDirection="row"
      border="1px solid gray"
    >
      {colorMap && <LegendGradient colorMap={colorMap} getValueLabel={getValueLabel} />}
    </Box>
    <Box height={10} position="relative">
      {colorMap && (
        <>
          <Box position="absolute" left={0}>
            <Typography>{getValueLabel(range[0])}</Typography>
          </Box>
          <Box position="absolute" right={0}>
            <Typography>{getValueLabel(range[1])}</Typography>
          </Box>
        </>
      )}
    </Box>
  </Box>
);
