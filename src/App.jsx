import React, { Fragment } from 'react';
import { Route, BrowserRouter as Router, Switch } from 'react-router-dom'

import Nav from './Nav'
import Map from './Map'

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
          <Map
            map_style="overview"
            dataSources={['road','bridges','rail', 'air', 'water']}
            dataLayers={[
              {key:'road', label: 'Roads'},
              {key:'bridges', label: 'National-roads bridges'},
              {key:'rail', label: 'Railways'},
              {key:'air', label: 'Airports'},
              {key:'water', label: 'Water'}
            ]}
            tooltipLayerSources={['road','bridges','rail', 'air', 'water']}
            />
        </Route>
        <Route path="/roads">
          <Map
            map_style="roads"
            dataSources={['road','bridges']}
            dataLayers={[
              {key:'road', label: 'Roads'},
              {key:'bridges', label: 'National-roads bridges'}
            ]}
            tooltipLayerSources={['road','bridges','flood']}
            />
        </Route>
        <Route path="/rail">
          <Map
            map_style="rail"
            dataSources={['rail']}
            dataLayers={[
              {key:'rail', label: 'Railways'}
            ]}
            tooltipLayerSources={['rail','flood']}
            />
        </Route>
        <Route path="/airwater">
          <Map
            map_style="airwater"
            dataSources={['air', 'water']}
            dataLayers={[
              {key:'air', label: 'Airports'},
              {key:'water', label: 'Water'}
            ]}
            tooltipLayerSources={['air','water','flood']}
            />
        </Route>
        <Route path="/flood">
          <Map
            map_style="flood"
            dataSources={[]}
            dataLayers={[]}
            tooltipLayerSources={['flood']}
            />
        </Route>
        <Route path="/adaptation">
          <Map
            zoom={6}
            lng={-61.5}
            lat={-34.6}
            map_style="adaptation"
            dataSources={['road','bridges']}
            dataLayers={[
              {key:'road', label: 'Roads'},
              {key:'bridges', label: 'National-roads bridges'}
            ]}
            tooltipLayerSources={['road','bridges']}
            />
        </Route>
        <Route path="/">
          <div className="jumbotron welcome-float">
            <h1 className="h1">Argentina Transport Risk Analysis - Results inquirer</h1>
            <p className="lead">

            This tool presents results of the modelling and analysis of climate-related risks
            to transport networks in Argentina.

            </p>
            <p className="lead">

            The modelling and analysis aim to support decision-making by identifying spatial
            criticailities, risks, and the performance of adaptation options under current and
            future fluvial and pluvial flooding outlooks. It comprises a network flow model,
            generation of failure scenarios, economic impact assessment, and cost-benefit
            analysis of adaptation options.

            </p>
            <p className="lead">

            The outputs visualized here were generated from a model created and documented
            here:

            </p>
            <p className="lead">

            <a href="https://github.com/oi-analytics/argentina-transport" target="blank">GitHub resources</a>

            </p>
            <p className="lead">

            <a href="https://argentina-transport-risk-analysis.readthedocs.io/en/latest/?badge=latest" target="blank">ReadTheDocs resources</a>

            </p>

            <h1 className="h1">Funding support</h1>
            <p className="lead">

            This results inquirer tool has been developed for the Government of Argentina with
            funding support from the World Bank Group and Global Facility for Disaster
            Reduction and Recovery (GFDRR).

            </p>
          </div>
        </Route>
      </Switch>
      </main>
    </Fragment>
  </Router>
)

export default App;
