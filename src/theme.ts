import { Theme, createTheme } from '@mui/material/styles';
import '@mui/styles';

declare module '@mui/styles/defaultTheme' {
  // eslint-disable-next-line @typescript-eslint/no-empty-interface
  interface DefaultTheme extends Theme {}
}

export const theme = createTheme({
  palette: {
    primary: {
      main: '#3f51b5',
    },
  },
  components: {
    MuiToolbar: {
      styleOverrides: {
        root: {
          background:"linear-gradient(180deg, rgba(0,128,0,1) 0%, rgba(0,128,0,1) 5%, rgba(9,9,9,1) 5%, rgba(9,9,9,1) 100%);"
        }
      }
    }
  }
});

export const globalStyleVariables = {
  sidebarWidth: 400,
  navbarHeight: 64,
  detailsSidebarWidth: 400,
};
