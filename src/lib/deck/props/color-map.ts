import { colorCssToRgb } from 'lib/helpers';
import _ from 'lodash';
import { Accessor, mergeTriggers, withTriggers } from './getters';

const memoizedColorCssToRgb = _.memoize(colorCssToRgb);

export function dataColorMap(dataSource: Accessor<any>, colorSource: Accessor<any>) {
  return withTriggers((x) => memoizedColorCssToRgb(colorSource(dataSource(x))), mergeTriggers(dataSource, colorSource));
}
