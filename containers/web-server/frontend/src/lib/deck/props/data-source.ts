import _ from 'lodash';
import { MapboxGeoJSONFeature } from 'mapbox-gl';

import { DataLoader } from '@/lib/data-loader/data-loader';

import { Accessor, withLoaderTriggers, withTriggers } from './getters';

export const featureProperty = _.memoize(
  (field: string | Accessor<any, MapboxGeoJSONFeature>): Accessor<any, MapboxGeoJSONFeature> => {
    return typeof field === 'string' ? withTriggers((f) => f.properties[field], [field]) : field;
  },
);

export function extraProperty(dataLoader: DataLoader): Accessor<any> {
  return withLoaderTriggers((f) => dataLoader.getData(f.id), dataLoader);
}
