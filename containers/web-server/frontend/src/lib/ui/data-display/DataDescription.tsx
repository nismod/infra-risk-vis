import { FC, useMemo } from 'react';

import { colorMap } from '@/lib/color-map';
import { ColorMap, ViewLayer } from '@/lib/data-map/view-layers';
import { DataItem } from '@/lib/ui/data-display/DataItem';

import { ColorBox } from './ColorBox';

export const DataDescription: FC<{
  viewLayer: ViewLayer;
  feature: any;
  colorMap: ColorMap;
}> = ({ viewLayer, feature, colorMap: { fieldSpec: colorField, colorSpec } }) => {
  const accessor = useMemo(() => viewLayer.dataAccessFn?.(colorField), [viewLayer, colorField]);

  const value = accessor?.(feature);

  const colorFn = useMemo(() => colorMap(colorSpec), [colorSpec]);

  const color = colorFn(value);

  const { getDataLabel, getValueFormatted } = viewLayer.dataFormatsFn(colorField);

  const dataLabel = getDataLabel(colorField);
  const formattedValue = getValueFormatted(value, colorField);

  return (
    <DataItem
      label={dataLabel}
      value={
        <>
          <ColorBox color={color} />
          {formattedValue ?? '-'}
        </>
      }
    />
  );
};
