import { MAP_SHAPE_TYPES, MapShapeType, shapeUrls } from './shapes';
import { load } from '@loaders.gl/core';
import { ImageLoader } from '@loaders.gl/images';

const iconHeight = 20,
  iconWidth = 20;

// one pixel more per each icon to avoid sampling pixels from neighboring icons
const iconWidthAtlas = 21;

const iconAtlas = (async function () {
  const n = MAP_SHAPE_TYPES.length;
  const canvasWidth = iconWidthAtlas * n;
  const images = await Promise.all(MAP_SHAPE_TYPES.map((x) => load(shapeUrls[x], ImageLoader, {})));

  const offscreen = new OffscreenCanvas(canvasWidth, iconHeight);
  const ctx = offscreen.getContext('2d');

  images.map((img, i) => ctx.drawImage(img, i * iconWidthAtlas, 0));

  return ctx.getImageData(0, 0, canvasWidth, iconHeight);
})();

const iconMapping = Object.fromEntries(
  MAP_SHAPE_TYPES.map((shapeType, i) => [
    shapeType,
    {
      x: i * iconWidthAtlas,
      y: 0,
      width: iconWidth,
      height: iconHeight,
      mask: true,
    },
  ]),
);

export function iconType(getIcon: MapShapeType | ((x: any) => MapShapeType)) {
  const iconGetter = typeof getIcon === 'function' ? (x: any) => getIcon(x) : () => getIcon;

  return {
    pointType: 'icon',
    iconAtlas,
    iconMapping,
    getIcon: iconGetter,
  };
}
