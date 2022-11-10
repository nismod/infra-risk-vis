import { Box } from '@mui/system';
import { FC, ReactNode, useMemo } from 'react';

import { ColorBox } from '@/lib/ui/data-display/ColorBox';
import { DataItem } from '@/lib/ui/data-display/DataItem';

import { formatRangeTruncation } from '../legend/GradientLegend';
import { ColorValue, RasterColorMapValues } from '../legend/RasterLegend';

export type RGBAColor = [number, number, number, number];

export function serializeColor([r, g, b, a]: [r: number, g: number, b: number, a: number]) {
  return `rgb(${r},${g},${b})`;
}
function useRasterColorMapLookup(colorMapValues: ColorValue[]): Record<string, { value: any; i: number }> {
  return useMemo(
    () => colorMapValues && Object.fromEntries(colorMapValues.map(({ value, color }, i) => [color, { value, i }])),
    [colorMapValues],
  );
}

export interface RasterBaseHoverProps {
  colorMap: RasterColorMapValues;
  color: RGBAColor;
  label: string;
  formatValue: (x: any) => ReactNode | string;
}

export const RasterBaseHover: FC<RasterBaseHoverProps> = ({ colorMap, color, label, formatValue }) => {
  const { colorMapValues, rangeTruncated = [false, false] } = colorMap;
  const rasterValueLookup = useRasterColorMapLookup(colorMapValues);

  const colorString = serializeColor(color);
  const valueLookup = rasterValueLookup?.[colorString];
  const { value, i } = valueLookup ?? {};
  return (
    <Box>
      <DataItem
        label={label}
        value={
          <>
            <ColorBox color={colorString} />
            {value == null ? formatValue(value) : formatRangeTruncation(formatValue(value), i, rangeTruncated)}
          </>
        }
      />
    </Box>
  );
};
