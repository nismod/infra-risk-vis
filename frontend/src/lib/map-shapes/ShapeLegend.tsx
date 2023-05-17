import { Box } from '@mui/material';
import { MapShapeType, shapeComponents } from './shapes';
export type { MapShapeType as LegendShapeType };

export const ShapeLegend = ({ type, color }: { type: MapShapeType; color: string }) => {
  const ShapeComponent = shapeComponents[type];
  return (
    <Box component="span" display="inline-block" marginRight={1}>
      <ShapeComponent fill={color} width="10px" height="10px" />
    </Box>
  );
};
