import { useMemo } from 'react';
import { useFetch } from 'use-http';

export function useRasterColorMapValues(colorScheme: string, stretchRange: [number, number]) {
  const [rangeMin, rangeMax] = stretchRange;

  const {
    loading,
    error,
    data: { colormap: colorMapValues = null } = {},
  } = useFetch(`/raster/colormap?colormap=${colorScheme}&stretch_range=[${rangeMin},${rangeMax}]`, { persist: true }, [
    colorScheme,
    rangeMin,
    rangeMax,
  ]);

  const result = useMemo(
    () => colorMapValues?.map(({ value, rgba: [r, g, b] }) => ({ value, color: `rgb(${r},${g},${b})` })),
    [colorMapValues],
  );

  return {
    loading,
    error,
    colorMapValues: result,
  };
}
