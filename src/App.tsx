import React from 'react';
import { Route, BrowserRouter as Router, Switch } from 'react-router-dom';
import CssBaseline from '@material-ui/core/CssBaseline';
import Toolbar from '@material-ui/core/Toolbar';

import Nav from './Nav';
import PageIntro from './PageIntro';
import { MapView } from './MapView';

import './index.css';
import 'mapbox-gl/dist/mapbox-gl.css';

const App = () => {
  return (
    <Router>
      <CssBaseline />
      <Nav />
      <Switch>
        <Route path="/" exact>
          <Toolbar /> {/* Prevents app bar from concealing content*/}
          <PageIntro />
        </Route>
        <Route path="/overview">
          <MapView />
        </Route>
      </Switch>
    </Router>
  );
};

export default App;
