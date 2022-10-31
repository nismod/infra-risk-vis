import { Box, Typography } from '@mui/material';
import { FC } from 'react';

const shapeStyles = {
  line: {
    borderRadius: 0,
    height: '6px',
    width: '12px',
  },
  circle: {
    width: '10px',
    height: '10px',
    borderRadius: '50%',
  },
  polygon: {
    width: '10px',
    height: '10px',
    borderRadius: '0',
  },
};

export type LayerLabelShapeType = keyof typeof shapeStyles;

export const LayerShapeLegend = ({ type, color }: { type: LayerLabelShapeType; color: string }) => {
  const shapeStyle = shapeStyles[type];
  return (
    <Box
      component="span"
      display="inline-block"
      marginRight="5px"
      marginBottom="1px"
      bgcolor={color}
      textAlign="center"
      {...shapeStyle}
    />
  );
};

export interface LayerLabelProps {
  label: string;
  type: LayerLabelShapeType;
  color: string;
}

export const LayerLabel: FC<LayerLabelProps> = ({ label, type, color }) => {
  return (
    <>
      <Typography>
        <LayerShapeLegend type={type} color={color} />
        {label}
      </Typography>
    </>
  );
};
