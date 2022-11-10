import { FC, ReactNode } from 'react';

import { useRasterColorMapValues } from '@/map/legend/use-color-map-values';

import { RasterColorMap } from '../legend/RasterLegend';
import { RasterBaseHover } from './RasterBaseHover';

export interface RasterHoverDescriptionProps {
  colorMap: RasterColorMap;
  color: [number, number, number, number];
  label: string;
  formatValue: (x: any) => ReactNode | string;
}

export const RasterHoverDescription: FC<RasterHoverDescriptionProps> = ({ colorMap, ...otherProps }) => {
  const { scheme, range, rangeTruncated } = colorMap;
  const colorMapValues = useRasterColorMapValues(scheme, range);
  const colorMapSpec = {
    colorMapValues,
    rangeTruncated,
  };

  return <RasterBaseHover colorMap={colorMapSpec} {...otherProps} />;
};
