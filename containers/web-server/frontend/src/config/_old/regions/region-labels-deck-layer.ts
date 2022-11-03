import { mvtLayer } from '@/lib/deck/layers/base';

import { BackgroundName } from '@/config/backgrounds';

import { REGIONS_METADATA, RegionLevel } from './metadata';

const LIGHT_TEXT = [240, 240, 240, 255];
const DARK_TEXT = [90, 90, 90, 255];

export function regionLabelsDeckLayer(level: RegionLevel, background: BackgroundName) {
  const config = REGIONS_METADATA[level];

  const color = background === 'satellite' ? LIGHT_TEXT : DARK_TEXT;
  return (
    config.showLabels &&
    mvtLayer({
      id: `boundaries_${level}-text`,
      data: `/vector/data/regions_${level}_labels.json`,
      loadOptions: {
        mvt: {
          layers: ['labels'],
        },
      },
      binary: false,
      minZoom: config.minZoom,
      pointType: 'text',
      getText: (f) => f.properties[config.fieldName],
      getTextSize: 24,
      getTextColor: color,
      textFontFamily: 'Arial',
      textFontWeight: 'bold',
      getPolygonOffset: ({ layerIndex }) => [0, -layerIndex * 100 - 2000],

      // won't work before deck.gl v8.7.0 is released (textFontSettings isn't mapped correctly)
      // see https://github.com/visgl/deck.gl/pull/6336
      //
      // textOutlineColor: [255, 255, 255, 255],
      // textOutlineWidth: 1,
      // textFontSettings: {
      //   sdf: true,
      // },
    })
  );
}
