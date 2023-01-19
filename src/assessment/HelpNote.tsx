import { FC } from 'react';
import { Paper, Typography } from '@mui/material';

export const HelpNote: FC<{ sx?: object }> = ({ sx, children }) => {
  const defaultSx = { mb: 1, p: 1, backgroundColor: '#fafafa', color: '#717171' };
  const actualSx = sx ? { ...defaultSx, ...sx } : defaultSx;
  return (
    <Paper sx={actualSx} elevation={0.5}>
      <Typography variant="body2">{children}</Typography>
    </Paper>
  );
};
