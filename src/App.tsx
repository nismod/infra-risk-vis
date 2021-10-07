import React from 'react';
import { Route, BrowserRouter as Router, Switch } from 'react-router-dom';
import { CssBaseline, Toolbar } from '@material-ui/core';

import { Nav } from './Nav';
import { IntroPage } from './IntroPage';
import { MapPage } from './MapPage';

import './index.css';
import 'mapbox-gl/dist/mapbox-gl.css';
import { ViewName } from './config/views';

export const App = () => {
  return (
    <Router>
      <CssBaseline />
      <Nav />
      <Switch>
        <Route path="/" exact>
          <Toolbar /> {/* Prevents app bar from concealing content*/}
          <IntroPage />
        </Route>
        <Route path="/:view(overview)" render={({ match: { params } }) => <MapPage view={params.view as ViewName} />} />
      </Switch>
    </Router>
  );
};
