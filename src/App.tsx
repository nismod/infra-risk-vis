import React from 'react';
import { Route, BrowserRouter as Router, Switch } from 'react-router-dom';
import { CssBaseline, StyledEngineProvider, Toolbar } from '@mui/material';
import { ThemeProvider } from '@mui/material/styles';

import { Nav } from './Nav';
import { IntroPage } from './IntroPage';
import { MapPage } from './MapPage';

import './index.css';
import 'mapbox-gl/dist/mapbox-gl.css';
import { ViewName } from './config/views';
import { RecoilRoot } from 'recoil';
import { theme } from './theme';

export const App = () => {
  return (
    <RecoilRoot>
      <StyledEngineProvider injectFirst>
        <ThemeProvider theme={theme}>
          <Router>
            <CssBaseline />
            <Nav />
            <Switch>
              <Route path="/" exact>
                <Toolbar /> {/* Prevents app bar from concealing content*/}
                <IntroPage />
              </Route>
              <Route
                path="/:view(overview)"
                render={({ match: { params } }) => <MapPage view={params.view as ViewName} />}
              />
            </Switch>
          </Router>
        </ThemeProvider>
      </StyledEngineProvider>
    </RecoilRoot>
  );
};
