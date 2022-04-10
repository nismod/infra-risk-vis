import { DataLoader } from 'lib/data-loader/data-loader';
import { MapboxGeoJSONFeature } from 'mapbox-gl';

export interface DataLoaderOptions {
  dataLoader: DataLoader;
}

// const mockDataRecord = {};

export function dataLoaderLayer(tileProps, { dataLoader }: DataLoaderOptions) {
  const {
    tile: { content },
  } = tileProps;
  if (content) {
    const ids: number[] = content.map((f: MapboxGeoJSONFeature) => f.id);

    dataLoader.loadDataForIds(ids);

    // mockDataRecord[z] ??= {};
    // mockDataRecord[z][x] ??= {};
    // mockDataRecord[z][x][y] ??= Object.fromEntries(content.map((f) => [f.properties.id, Math.random() * 1_000_000]));
  }

  return null;
}
