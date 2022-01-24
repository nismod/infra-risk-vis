import { MVTLayer } from 'deck.gl';

export function infrastructureLayer(...props) {
  return new MVTLayer(
    {
      binary: true,
      autoHighlight: true,
      highlightColor: [0, 255, 255, 255],
      refinementStrategy: 'best-available',
    } as any,
    ...props,
  );
}
