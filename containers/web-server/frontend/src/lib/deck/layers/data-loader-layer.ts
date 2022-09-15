import { DataLoader } from 'lib/data-loader/data-loader';
import { MapboxGeoJSONFeature } from 'mapbox-gl';

export interface DataLoaderOptions {
  dataLoader: DataLoader;
}

export function dataLoaderLayer(tileProps, { dataLoader }: DataLoaderOptions) {
  const {
    tile: { content },
  } = tileProps;
  if (content && dataLoader) {
    const ids: number[] = content.map((f: MapboxGeoJSONFeature) => f.id);

    dataLoader.loadDataForIds(ids);
  }

  return null;
}
