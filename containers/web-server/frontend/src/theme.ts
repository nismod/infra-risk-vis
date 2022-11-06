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
    MuiCssBaseline: {
      styleOverrides: `
        @font-face {
          font-family: 'Catamaran';
          src: url('/fonts/static/catamaran-regular.ttf');
          font-weight: normal;
        }
        @font-face {
          font-family: 'Catamaran';
          src: url('/fonts/static/catamaran-semibold.ttf');
          font-weight: 600;
        }
        @font-face {
          font-family: 'Catamaran';
          src: url('/fonts/static/catamaran-bold.ttf');
          font-weight: bold;
        }
      `,
    },
    MuiAccordion: {
      styleOverrides: {
        root: {
          /**
           * Change disabled accordion style to white background,
           * by default text will still be grayed out
           */
          '&.Mui-disabled': {
            backgroundColor: 'white',
          },
        },
      },
    },
    MuiTableRow: {
      styleOverrides: {
        root: {
          verticalAlign: 'top',
        }
      }
    }
  },
  typography: {
    h1: {
      fontWeight: 600,
      fontFamily: 'Catamaran, Roboto, sans-serif',
      fontSize: '2.75rem',
      margin: '1rem 0',
    },
    h2: {
      fontWeight: 600,
      fontFamily: 'Catamaran, Roboto, sans-serif',
      fontSize: '1.6rem',
      margin: '0.5rem 0',
    },
    h3: {
      fontWeight: 600,
      fontFamily: 'Catamaran, Roboto, sans-serif',
      fontSize: '1.2rem',
    },
    h6: {
      fontWeight: 600,
      fontFamily: 'Catamaran, Roboto, sans-serif',
      fontSize: '1.6rem',
    },
  },
});

export const globalStyleVariables = {
  controlSidebarWidth: 400,
  detailSidebarWidth: 500,
  navbarHeight: 64,
  detailsSidebarWidth: 400,
};
