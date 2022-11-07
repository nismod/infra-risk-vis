import { Box } from '@mui/material';
import { FC, useMemo } from 'react';

import { ColorBox } from '@/lib/ui/data-display/ColorBox';
import { DataItem } from '@/lib/ui/data-display/DataItem';

import { useRasterColorMapValues } from '@/map/legend/use-color-map-values';

function useRasterColorMapLookup(colorMapValues) {
  return useMemo(
    () => colorMapValues && Object.fromEntries(colorMapValues.map(({ value, color }) => [color, value])),
    [colorMapValues],
  );
}

export interface RasterHoverDescriptionProps {
  scheme: string;
  range: [number, number];
  color: [number, number, number, number];
  label: string;
  formatValue: (x: any) => string;
}

export const RasterHoverDescription: FC<RasterHoverDescriptionProps> = ({
  scheme,
  range,
  color,
  label,
  formatValue,
}) => {
  const colorMapValues = useRasterColorMapValues(scheme, range);
  const rasterValueLookup = useRasterColorMapLookup(colorMapValues);

  const colorString = `rgb(${color[0]},${color[1]},${color[2]})`;
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
