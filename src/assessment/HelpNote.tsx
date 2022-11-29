import { Paper, Typography } from '@mui/material';

export const HelpNote = ({ children }) => {
  return (
    <Paper sx={{ mb: 1, p: 1, backgroundColor: '#fafafa', color: '#717171' }} elevation={0.5}>
      <Typography variant="body2">{children}</Typography>
    </Paper>
  );
};
