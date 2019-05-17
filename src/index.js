import React, { Fragment } from 'react';
import { render } from 'react-dom'
import { Route, BrowserRouter as Router, Switch } from 'react-router-dom'

import Nav from './Nav'
import StaticMap from './StaticMap'
import TooltipMap from './TooltipMap'

import './index.css';
import 'bootstrap/dist/css/bootstrap.min.css'

const App = () => (
  <Router>
    <Fragment>
      <Route path="/" component={Nav}/>
      <main class="map-height">
      <Switch>
        <Route path="/overview">
          <StaticMap
            map_style={"http://localhost:8080/styles/overview/style.json"}
            toggleableLayerIds={['road_secondary_tertiary', 'road_trunk_primary', 'road_major_motorway', 'aeroway', 'waterway']}
            clickableLayerAttributes={{
              'road_secondary_tertiary': {
                '_header': 'Rural Road',
                'gml_id': 'Name',
                'u_jurisdic': 'Classification'
              },
              'road_trunk_primary': {
                '_header': 'Provincial Road',
                'name': 'Name',
                'jurisdicti': 'Classification'
              },
              'road_major_motorway': {
                '_header': 'Motorway',
                'name': 'Name',
                'jurisdicti': 'Classification'
              },
              'waterway': {},
              'aeroway': {
                '_header': 'Flight path',
                'from_iata': 'From',
                'to_iata': 'To'
              }
            }}/>
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
