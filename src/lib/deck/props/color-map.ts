import { colorCssToRgb } from 'lib/helpers';
import { Accessor, mergeTriggers, withTriggers } from './getters';

export function dataColorMap(dataSource: Accessor<any>, colorSource: Accessor<any>) {
  return withTriggers(
    (x) => colorCssToRgb(colorSource(dataSource(x))),
    mergeTriggers(dataSource, colorSource),
  );
}
