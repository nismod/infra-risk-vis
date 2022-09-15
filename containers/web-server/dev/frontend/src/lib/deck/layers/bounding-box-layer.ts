import { PathStyleExtension } from '@deck.gl/extensions';

import { geoJsonLayer } from './base';

export interface BoundingBoxLayerOptions {
  bboxGeom: any;
}

export function boundingBoxLayer({ bboxGeom }: BoundingBoxLayerOptions, ...props) {
  return geoJsonLayer(
    {
      data: bboxGeom,

      stroked: true,
      filled: false,
      getLineColor: [0, 255, 255],
      lineWidthMinPixels: 1,
      getLineWidth: 2,
      lineWidthUnits: 'pixels',

      getDashArray: [5, 2],
      dashJustified: true,
      dashGapPickable: true,
      extensions: [new PathStyleExtension({ dash: true })],
    },
    props,
  );
}
