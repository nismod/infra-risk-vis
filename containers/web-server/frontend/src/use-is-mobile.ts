import { useMediaQuery } from '@mui/material';

export function useIsMobile() {
  return useMediaQuery((theme: any) => theme.breakpoints.down('md'));
}
