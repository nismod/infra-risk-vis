import { Button } from '@mui/material';

import { withProps } from 'lib/react/with-props';

declare module '@mui/material/Button' {
  interface ButtonPropsColorOverrides {
    map: true;
  }
}

export const MapHudButton = withProps(
  Button,
  {
    variant: 'contained',
    color: 'map',
    sx: (theme) => ({
      paddingInline: 0,
      minWidth: 'auto',
      width: '40px',
      height: '36px',
      '&.Mui-disabled': {
        opacity: 1,
        color: theme.palette.text.disabled,
      },
    }),
  },
  'MapHudButton',
);
