import bbox from '@turf/bbox';
import bboxPolygon from '@turf/bbox-polygon';
import buffer from '@turf/buffer';

/**
 * 2D bbox format as defined in GeoJSON, turf etc:
 * [minX, minY, maxX, maxY]
 */
export type BoundingBox = [
  number, // minX
  number, // minY
  number, // maxX
  number, // maxY
];

/**
 * Nominatim search result bbox format:
 * [minY, maxY, minX, maxX]
 */
export type NominatimBoundingBox = [
  number, // minY
  number, // maxY
  number, // minX
  number, // maxX
];

/**
 * Deck.GL bbox format:
 * [[minX, minY], [maxX, maxY]]
 */
export type DeckBoundingBox = [
  [
    number, // minX
    number, // minY
  ],
  [
    number, // maxX
    number, // maxY
  ],
];

export function appToDeckBoundingBox(appBoundingBox: BoundingBox): DeckBoundingBox {
  return [
    [appBoundingBox[0], appBoundingBox[1]],
    [appBoundingBox[2], appBoundingBox[3]],
  ];
}

export function nominatimToAppBoundingBox(nominatimBoundingBox: NominatimBoundingBox): BoundingBox {
  return [nominatimBoundingBox[2], nominatimBoundingBox[0], nominatimBoundingBox[3], nominatimBoundingBox[1]];
}

export function extendBbox(boundingBox: BoundingBox, kilometers: number): BoundingBox {
  return bbox(buffer(bboxPolygon(boundingBox), kilometers)) as BoundingBox;
}
