import React from 'react';
import { Route, BrowserRouter as Router, Switch } from 'react-router-dom';
import CssBaseline from '@material-ui/core/CssBaseline';
import Toolbar from '@material-ui/core/Toolbar';

import Nav from './Nav';
import Map from './Map';
import PageIntro from './PageIntro';

import './index.css';
import 'mapbox-gl/dist/mapbox-gl.css';

const App = () => (
  <Router>
    <CssBaseline />
    <Route path="/" component={Nav}/>
    <Switch>
      <Route path="/" exact>
        <Toolbar />
        <PageIntro />
      </Route>
      <Route path="/overview">
        <Map
          map_style="overview"
          dataSources={[
            'electricity',
            'rail',
            'roads_main',
            'roads_other'
          ]}
          dataLayers={[
            {key: 'electricity', label: 'Power Grid', linear: true, color: "#eca926"},
            {key: 'rail', label: 'Railways', linear: true, color: "#444"},
            {key: 'trunk', label: 'Trunk Roads', linear: true, color: "#941339"},
            {key: 'motorway', label: 'Motorways', linear: true, color: "#941339"},
            {key: 'primary', label: 'Primary Roads', linear: true, color: "#cb3e4e"},
            {key: 'secondary', label: 'Secondary Roads', linear: true, color: "#8471a8"},
            {key: 'roads_other', label: 'Tertiary and Other Roads', linear: true, color: "#b2afaa"},

          ]}
          tooltipLayerSources={[
            'electricity',
            'rail',
            'roads_main',
            'roads_other'
          ]}
          />
      </Route>
    </Switch>
  </Router>
)

export default App;
