import { FC } from 'react';

import { useRasterColorMapValues } from '@/map/legend/use-color-map-values';

import { RasterBaseHover } from './RasterBaseHover';

export interface RasterHoverDescriptionProps {
  scheme: string;
  range: [number, number];
  color: [number, number, number, number];
  label: string;
  formatValue: (x: any) => string;
}

export const RasterHoverDescription: FC<RasterHoverDescriptionProps> = ({ scheme, range, ...otherProps }) => {
  const colorMapValues = useRasterColorMapValues(scheme, range);

  return <RasterBaseHover colorMapValues={colorMapValues} {...otherProps} />;
};
