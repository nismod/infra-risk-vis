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
        <Route path="/roads">
          <StaticMap
            map_style={"http://localhost:8080/styles/roads/style.json"}
            dataSources={['road','bridges']}
            dataLayers={[
              {key:'road', label: 'Roads'},
              {key:'bridges', label: 'National-roads bridges'}
            ]}
            />
        </Route>
        <Route path="/rail">
          <StaticMap
            map_style={"http://localhost:8080/styles/rail/style.json"}
            dataSources={['rail']}
            dataLayers={[
              {key:'rail', label: 'Railways'}
            ]}
            />
        </Route>
        <Route path="/airwater">
          <StaticMap
            map_style={"http://localhost:8080/styles/airwater/style.json"}
            dataSources={['air', 'water']}
            dataLayers={[
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
            <h1 className="h1">Argentina Transport Risk Analysis - Results inquirer</h1>
            <p className="lead">
            This tool presents results of the modelling and analysis of climate-related risks to transport networks in Argentina.
            </p>
            <p className="lead">
            The modelling and analysis aim to support decision-making by identifying spatial criticailities, risks, and the performance 
            of adaptation options under current and future fluvial and pluvial flooding outlooks. It comprises a network flow model, 
            generation of failure scenarios, economic impact assessment, and cost-benefit analysis of adaptation options.
            </p>
            <p className="lead">
            The outputs visalized here were generated from a model created and documented here:
            </p>
            <p className="lead">
            <a href="https://github.com/oi-analytics/argentina-transport" target="blank">Github resources</a>
            </p>
            <p className="lead">
            <a href="https://argentina-transport-risk-analysis.readthedocs.io/en/latest/?badge=latest" target="blank">Readthedocs resources</a>  
            </p>
            <h1 className="h1">Funding support</h1>
            <p className="lead">
            This results inquirer tool has been developed for the Government of Argentina with funding 
            support from the World Bank Group and Global Facility for Disaster Reduction and Recovery (GFDRR).   
            </p>
          </div>
        </Route>
      </Switch>
      </main>
    </Fragment>
  </Router>
)

render(<App />, document.getElementById('root'));
