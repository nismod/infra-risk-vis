import React, { Fragment } from 'react';
import { render } from 'react-dom'
import { Route, BrowserRouter as Router, Switch } from 'react-router-dom'

import Nav from './Nav'
import StaticMap from './StaticMap'
import TooltipMap from './TooltipMap'

import './index.css';
import 'bootstrap/dist/css/bootstrap.min.css'
import 'mapbox-gl/dist/mapbox-gl.css'

const App = () => (
  <Router>
    <Fragment>
      <Route path="/" component={Nav}/>
      <main className="map-height">
      <Switch>
        <Route path="/overview">
          <StaticMap
            map_style={"http://localhost:8080/styles/overview/style.json"}
            dataSources={['road', 'air', 'water']}
            dataLayers={['road_all', 'road_national', 'air', 'water']}
            />
          </Route>
          <Route path="/flood">
            <TooltipMap
              map_style={"http://localhost:8080/styles/flood/style.json"}
              tooltipLayerSources={['flood']}/>
          </Route>
      </Switch>
      </main>
    </Fragment>
  </Router>
)

render(<App />, document.getElementById('root'));
