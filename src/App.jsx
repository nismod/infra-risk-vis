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
            'road_edges',
            'bridges',
            'rail_edges',
            'rail_nodes',
            'elec_edges',
            'elec_nodes'
          ]}
          dataLayers={[
            {key: 'road_edges', label: 'Roads', linear: true, color: "#b2afaa"},
            {key: 'bridges', label: 'Bridges', color: "#487dbc"},
            {key: 'rail_edges', label: 'Railways', linear: true, color: "#444"},
            {key: 'rail_nodes', label: 'Stations', color: "#444"},
            {key: 'elec_edges_high', label: 'Power Lines (High Voltage)', linear: true, color: "#eca926"},
            {key: 'elec_edges_low', label: 'Power Lines (Low Voltage)', linear: true, color: "#f1d75c"},
            {key: 'elec_nodes', label: 'Power Nodes', color: "#eca926"},
          ]}
          tooltipLayerSources={[
            'road_edges',
            'bridges',
            'rail_edges',
            'rail_nodes',
            'elec_edges',
            'elec_nodes'
          ]}
          />
      </Route>
    </Switch>
  </Router>
)

export default App;
