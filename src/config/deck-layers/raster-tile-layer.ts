import { TileLayer, BitmapLayer } from 'deck.gl';

function getBoundsForTile(tileProps) {
  const {
    bbox: { west, south, east, north },
  } = tileProps;

  return [west, south, east, north];
}

export function rasterTileLayer(bitmapProps, ...props) {
  return new TileLayer(...props, {
    renderSubLayers: (tileProps) =>
      new BitmapLayer(
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
