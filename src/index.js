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
            dataSources={['road','bridges','rail', 'air', 'water']}
            dataLayers={[
              {key:'road', label: 'Roads'},
              {key:'bridges', label: 'National-roads bridges'},
              {key:'rail', label: 'Railways'},
              {key:'air', label: 'Airports'}, 
              {key:'water', label: 'Water'}
            ]}
            />
        </Route>
        <Route path="/flood">
          <TooltipMap
            map_style={"http://localhost:8080/styles/flood/style.json"}
            tooltipLayerSources={['flood']}/>
        </Route>
        <Route path="/">
          <div className="jumbotron welcome-float">
            <h1 className="h1">Argentina Transport Risk Analysis</h1>
            <p className="lead">
            This tool presents results of the modelling and analysis of climate-related risks to transport networks in Argentina.
            </p>
            <p className="lead">
            The modelling and analysis aim to support decision-making by identifying the performance of adaptation options under current and future scenarios. It comprises a network flow model, generation of failure scenarios, economic impact assessment, and cost-benefit analysis of adaptation options.
            </p>
          </div>
        </Route>
      </Switch>
      </main>
    </Fragment>
  </Router>
)

render(<App />, document.getElementById('root'));
