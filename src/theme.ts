import { Theme, createTheme } from '@mui/material/styles';
import '@mui/styles';

declare module '@mui/styles/defaultTheme' {
  // eslint-disable-next-line @typescript-eslint/no-empty-interface
  interface DefaultTheme extends Theme {}
}

export const theme = createTheme({
  palette: {
    text: {
      primary: '#222',
    },
    primary: {
      main: '#172617',
    },
  },
  components: {
    MuiCssBaseline: {
      styleOverrides: `
      @font-face {
        font-family: 'Inter';
        font-style:normal;
        font-weight: 400;
        font-display:swap;
        src: url('/fonts/Inter-Regular.woff2') format("woff2");
      }
      @font-face {
        font-family: 'InterDisplay';
        font-style:normal;
        font-weight: 400;
        font-display:swap;
        src: url('/fonts/InterDisplay-Regular.woff2') format("woff2");
      }

      @font-face {
        font-family: 'Inter';
        font-style:normal;
        font-weight: 500;
        font-display:swap;
        src: url('/fonts/Inter-Medium.woff2') format("woff2");
      }
      @font-face {
        font-family: 'InterDisplay';
        font-style:normal;
        font-weight: 500;
        font-display:swap;
        src: url('/fonts/InterDisplay-Medium.woff2') format("woff2");
      }

      @font-face {
        font-family: 'Inter';
        font-style:normal;
        font-weight: 600;
        font-display:swap;
        src: url('/fonts/Inter-SemiBold.woff2') format("woff2");
      }
      @font-face {
        font-family: 'InterDisplay';
        font-style:normal;
        font-weight: 600;
        font-display:swap;
        src: url('/fonts/InterDisplay-SemiBold.woff2') format("woff2");
      }
      `,
    },
  },
  typography: {
    h1: {
      fontWeight: 600,
      fontFamily: 'InterDisplay, sans-serif',
      fontSize: '3rem',
      '@media (min-width:900px)': {
        fontSize: '4.5rem',
      },
      margin: '1rem 0',
      letterSpacing: '-2px',
      lineHeight: 1.1,
      maxWidth: '11em',
    },
    h2: {
      fontWeight: 600,
      fontFamily: 'Inter, sans-serif',
      fontSize: '1.5rem',
      margin: '0.5rem 0',
    },
    h3: {
      fontWeight: 600,
      fontFamily: 'Inter, sans-serif',
      fontSize: '1.25rem',
    },
    h5: {
      fontWeight: 600,
      fontFamily: 'Inter, sans-serif',
      fontSize: '1.5rem',
      letterSpacing: '-1px',
      margin: '1rem 0 2rem',
      lineHeight: 1.2,
    },
    h6: {
      fontWeight: 600,
      fontFamily: 'Inter, sans-serif',
      fontSize: '1.5rem',
      letterSpacing: '-1px',
    },
  },
});

export const globalStyleVariables = {
  controlSidebarWidth: 400,
  detailSidebarWidth: 500,
  navbarHeight: 64,
  detailsSidebarWidth: 400,
};
