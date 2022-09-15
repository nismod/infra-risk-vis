import { useMemo } from 'react';
import { useFetch } from 'use-http';

export function useRasterColorMapValues(colorScheme: string, stretchRange: [number, number]) {
  const [rangeMin, rangeMax] = stretchRange;

  const {
    loading,
    error,
    data: { colormap: colorMapValues = null } = {},
  } = useFetch(`/raster/colormap?colormap=${colorScheme}&stretch_range=[0,1]`, { persist: false }, [colorScheme]);

  const rangeSize = rangeMax - rangeMin;

  const result = useMemo(
    () =>
      colorMapValues?.map(({ value, rgba: [r, g, b] }) => ({
        value: rangeMin + value * rangeSize,
        color: `rgb(${r},${g},${b})`,
      })),
    [colorMapValues, rangeMin, rangeSize],
  );

  return {
    loading,
    error,
    colorMapValues: result,
  };
}
