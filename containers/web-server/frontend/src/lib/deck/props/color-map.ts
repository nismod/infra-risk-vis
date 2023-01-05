import _ from 'lodash';

import { colorCssToRgb } from '@/lib/helpers';

import { Accessor, mergeTriggers, withTriggers } from './getters';

const memoizedColorCssToRgb = _.memoize(colorCssToRgb);

export function dataColorMap<T>(dataSource: Accessor<T>, colorSource: Accessor<string, T>) {
  return withTriggers(
    (x) => memoizedColorCssToRgb(colorSource(dataSource(x))),
    mergeTriggers(dataSource, colorSource),
  );
}
