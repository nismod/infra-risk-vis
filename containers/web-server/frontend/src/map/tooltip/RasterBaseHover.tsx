import { Box } from '@mui/system';
import { FC, useMemo } from 'react';

import { ColorBox } from '@/lib/ui/data-display/ColorBox';
import { DataItem } from '@/lib/ui/data-display/DataItem';

import { ColorValue } from '../legend/GradientLegend';

export type RGBAColor = [number, number, number, number];

export function serializeColor([r, g, b, a]: [r: number, g: number, b: number, a: number]) {
  return `rgb(${r},${g},${b})`;
}
function useRasterColorMapLookup(colorMapValues: ColorValue[]) {
  return useMemo(
    () => colorMapValues && Object.fromEntries(colorMapValues.map(({ value, color }) => [color, value])),
    [colorMapValues],
  );
}

export interface RasterBaseHoverProps {
  colorMapValues: ColorValue[];
  color: RGBAColor;
  label: string;
  formatValue: (x: any) => string;
}

export const RasterBaseHover: FC<RasterBaseHoverProps> = ({ colorMapValues, color, label, formatValue }) => {
  const rasterValueLookup = useRasterColorMapLookup(colorMapValues);

  const colorString = serializeColor(color);
  const value = rasterValueLookup?.[colorString];
  return (
    <Box>
      <DataItem
        label={label}
        value={
          <>
            <ColorBox color={colorString} />
            {formatValue(value)}
          </>
        }
      />
    </Box>
  );
};
