import { selectorFamily, useRecoilValue } from 'recoil';

import { apiClient } from '@/api-client';

import { ColorValue } from './GradientLegend';

const colorMapValuesQuery = selectorFamily({
  key: 'colorMapValuesQuery',
  get: (colorScheme: string) => async () => {
    return await apiClient.colormap.colormapGetColormap({
      colormap: colorScheme,
      stretchRange: '[0,1]',
    });
  },
});

const colorMapValuesState = selectorFamily({
  key: 'colorMapValuesState',
  get:
    (colorSpec: { scheme: string; range: [number, number] }) =>
    ({ get }) => {
      const values = get(colorMapValuesQuery(colorSpec.scheme));
      const [rangeMin, rangeMax] = colorSpec.range;

      const rangeSize = rangeMax - rangeMin;

      return values.colormap.map(({ value, rgba: [r, g, b] }) => ({
        value: rangeMin + value * rangeSize,
        color: `rgb(${r},${g},${b})`,
      }));
    },
});

export function useRasterColorMapValues(colorScheme: string, stretchRange: [number, number]): ColorValue[] {
  return useRecoilValue(colorMapValuesState({ scheme: colorScheme, range: stretchRange }));
}
