import { Box, Typography } from '@mui/material';

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

export const LayerShapeLegend = ({ type, color }: { type: 'line' | 'circle' | 'square'; color: string }) => {
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

export const LayerLabel = ({ label, type, color }) => {
  return (
    <>
      <Typography>
        <LayerShapeLegend type={type} color={color} />
        {label}
      </Typography>
    </>
  );
};
