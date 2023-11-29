import React from 'react';
import { RecoilRoot } from 'recoil';
import { Route, BrowserRouter as Router, Switch } from 'react-router-dom';
import { CssBaseline, StyledEngineProvider, Box } from '@mui/material';
import { ThemeProvider } from '@mui/material/styles';

import { Nav } from './Nav';
import { IntroPage } from './pages/IntroPage';
import { MapPage } from './pages/MapPage';
import { AssessmentPage } from 'pages/AssessmentPage';
import { DataPage } from './pages/DataPage';
import { theme, globalStyleVariables } from './theme';

import './index.css';
import 'mapbox-gl/dist/mapbox-gl.css';

export const App = () => {
  return (
    <RecoilRoot>
      <StyledEngineProvider injectFirst>
        <ThemeProvider theme={theme}>
          <Router>
            <CssBaseline />
            <Nav height={globalStyleVariables.navbarHeight} />
            <Box position="absolute" top={globalStyleVariables.navbarHeight} bottom={0} left={0} right={0}>
              <Switch>
                <Route path="/" exact>
                  <IntroPage />
                </Route>
                <Route
                  path="/:view(exposure|risk|adaptation|nature-based-solutions)"
                  render={({ match: { params } }) => <MapPage view={params.view} />}
                />
                <Route path="/assessment" exact>
                  <AssessmentPage />
                </Route>
                <Route path="/data" exact>
                  <DataPage />
                </Route>
              </Switch>
            </Box>
          </Router>
        </ThemeProvider>
      </StyledEngineProvider>
    </RecoilRoot>
  );
};
