import { Box, Typography } from '@mui/material';
import { FC, ReactNode } from 'react';

const legendHeight = 10;

export interface ColorValue {
  color: string;
  value: any;
}

const LegendGradient: FC<{
  colorMapValues: ColorValue[];
  getValueLabel: (value: number) => string;
}> = ({ colorMapValues, getValueLabel }) => {
  return (
    <>
      {colorMapValues.map(({ color, value }, i) => (
        <Box key={i} height={legendHeight} width={1} bgcolor={color} title={getValueLabel(value)} />
      ))}
    </>
  );
};

export interface GradientLegendProps {
  label: string | ReactNode;
  description?: string;
  range: [number, number];
  colorMapValues: ColorValue[];
  getValueLabel: (x: any) => string;
}

export const GradientLegend: FC<GradientLegendProps> = ({
  label,
  description,
  range,
  colorMapValues,
  getValueLabel,
}) => (
  <Box mb={2}>
    <Box mb={1}>
      <Typography variant="body1">{label}</Typography>
      {description && <Typography variant="body2">{description}</Typography>}
    </Box>
    <Box
      height={legendHeight + 2}
      width={255}
      bgcolor="#ccc"
      display="flex"
      flexDirection="row"
      border="1px solid gray"
    >
      {colorMapValues && <LegendGradient colorMapValues={colorMapValues} getValueLabel={getValueLabel} />}
    </Box>
    <Box height={10} position="relative">
      {colorMapValues && (
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
