import { bitmapLayer, tileLayer } from './base';

function getBoundsForTile(tileProps) {
  const {
    bbox: { west, south, east, north },
  } = tileProps;

  return [west, south, east, north];
}

export function rasterTileLayer(bitmapProps, ...props) {
  return tileLayer(props, {
    renderSubLayers: (tileProps) =>
      bitmapLayer(
        tileProps,
        {
          data: null,
          image: tileProps.data,
          bounds: getBoundsForTile(tileProps.tile),
        },
        bitmapProps,
      ),
  });
}
