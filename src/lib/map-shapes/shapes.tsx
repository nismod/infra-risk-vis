import { default as CircleShapeSrc, ReactComponent as CircleShape } from './shapes/circle.svg';
import { default as SquareShapeSrc, ReactComponent as SquareShape } from './shapes/square.svg';
import { default as PolygonShapeSrc, ReactComponent as PolygonShape } from './shapes/polygon.svg';
import { default as LineShapeSrc, ReactComponent as LineShape } from './shapes/line.svg';
import { default as InvTriangleShapeSrc, ReactComponent as InvTriangleShape } from './shapes/inv-triangle.svg';
import { default as DiamondShapeSrc, ReactComponent as DiamondShape } from './shapes/diamond.svg';
import _ from 'lodash';

export const MAP_SHAPE_TYPES = ['line', 'circle', 'square', 'polygon', 'inv-triangle', 'diamond'] as const;
export type MapShapeType = typeof MAP_SHAPE_TYPES[number];

type SVGComponent = typeof LineShape;

export const shapeComponents: Record<MapShapeType, SVGComponent> = {
  line: LineShape,
  circle: CircleShape,
  square: SquareShape,
  polygon: PolygonShape,
  'inv-triangle': InvTriangleShape,
  diamond: DiamondShape,
};

export const shapeUrls: Record<MapShapeType, string> = {
  line: LineShapeSrc,
  circle: CircleShapeSrc,
  square: SquareShapeSrc,
  polygon: PolygonShapeSrc,
  'inv-triangle': InvTriangleShapeSrc,
  diamond: DiamondShapeSrc,
};
