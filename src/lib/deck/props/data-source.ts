import _ from 'lodash';
import { MapboxGeoJSONFeature } from 'mapbox-gl';
import { Accessor, withTriggers } from './getters';

export const featureProperty = _.memoize(function (
  field: string | Accessor<any, MapboxGeoJSONFeature>,
): Accessor<any, MapboxGeoJSONFeature> {
  return typeof field === 'string'
    ? withTriggers((f) => f.properties[field], [field])
    : field;
});
