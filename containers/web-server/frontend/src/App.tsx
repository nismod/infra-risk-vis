import { Box, CssBaseline, StyledEngineProvider } from '@mui/material';
import { ThemeProvider } from '@mui/material/styles';
import { Route, BrowserRouter as Router, Switch } from 'react-router-dom';
import { RecoilRoot } from 'recoil';

import { Nav } from './Nav';
import { DataPage } from './pages/DataPage';
import { IntroPage } from './pages/IntroPage';
import { MapPage } from './pages/MapPage';
import { TermsPage } from './pages/TermsPage';
import { globalStyleVariables, theme } from './theme';

import 'mapbox-gl/dist/mapbox-gl.css';
import './index.css';

export const App = () => {
  return (
    <RecoilRoot>
      <StyledEngineProvider injectFirst>
        <ThemeProvider theme={theme}>
          <Router>
            <CssBaseline />
            <Nav height={globalStyleVariables.navbarHeight} />
            <Box
              position="absolute"
              top={globalStyleVariables.navbarHeight}
              bottom={0}
              left={0}
              right={0}
            >
              <Switch>
                <Route path="/" exact>
                  <IntroPage />
                </Route>
                <Route
                  path="/view/:view"
                  render={({ match: { params } }) => <MapPage view={params.view} />}
                />
                <Route path="/data" exact>
                  <DataPage />
                </Route>
                <Route path="/terms-of-use" exact>
                  <TermsPage />
                </Route>
              </Switch>
            </Box>
          </Router>
        </ThemeProvider>
      </StyledEngineProvider>
    </RecoilRoot>
  );
};
