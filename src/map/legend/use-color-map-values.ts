import { useFetch } from 'use-http';

export function useColorMapValues(colorScheme: string, stretchRange: [number, number]) {
  const [rangeMin, rangeMax] = stretchRange;

  const {
    loading,
    error,
    data: { colormap: colorMapValues = null } = {},
  } = useFetch(`http://localhost:5000/colormap?colormap=${colorScheme}&stretch_range=[${rangeMin},${rangeMax}]`, {}, [
    colorScheme,
    rangeMin,
    rangeMax,
  ]);

  return { loading, error, colorMapValues };
}
