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
            dataSources={[
              'road'
            ]}
            dataLayers={[
              {key:'road_class_1', label: 'Road Class 1', linear: true, color: "#000004"},
              {key:'road_class_2', label: 'Road Class 2', linear: true, color: "#2c115f"},
              {key:'road_class_3', label: 'Road Class 3', linear: true, color: "#721f81"},
              {key:'road_class_4', label: 'Road Class 4', linear: true, color: "#b73779"},
              {key:'road_class_5', label: 'Road Class 5', linear: true, color: "#f1605d"},
              {key:'road_class_6', label: 'Road Class 6', linear: true, color: "#feb078"},
              {key:'Cambodia_Energy', label: 'Cambodia Energy Network', linear: true, color: "#000000"},
              {key:'Laos_Energy', label: 'Lao PDR Energy Network', linear: true, color: "#000000"},
              {key:'Myanmar_Energy', label: 'Myanmar Energy Network', linear: true, color: "#000000"},
              {key:'Thailand_Energy', label: 'Thailand Energy Network', linear: true, color: "#000000"},
              {key:'Vietnam_Energy', label: 'Vietnam Energy Network', linear: true, color: "#000000"}
            ]}
            tooltipLayerSources={[
              'road'
            ]}
            />
        </Route>
        <Route path="/roads">
          <Map
            map_style="roads"
            dataSources={[
              'road'
            ]}
            dataLayers={[
              {key:'road_class_1', label: 'Road Class 1', linear: true, color: "#000004"},
              {key:'road_class_2', label: 'Road Class 2', linear: true, color: "#2c115f"},
              {key:'road_class_3', label: 'Road Class 3', linear: true, color: "#721f81"},
              {key:'road_class_4', label: 'Road Class 4', linear: true, color: "#b73779"},
              {key:'road_class_5', label: 'Road Class 5', linear: true, color: "#f1605d"},
              {key:'road_class_6', label: 'Road Class 6', linear: true, color: "#feb078"}

            ]}
            tooltipLayerSources={[
              'road'
            ]}
            />
        </Route>
        <Route path="/energy_network">
          <Map
            map_style="energy_network"
            dataSources={[
              'energy_network'
            ]}
            dataLayers={[
              {key:'Cambodia_Energy', label: 'Cambodia Energy Network', linear: true, color: "#000000"},
              {key:'Laos_Energy', label: 'Lao PDR Energy Network', linear: true, color: "#000000"},
              {key:'Myanmar_Energy', label: 'Myanmar Energy Network', linear: true, color: "#000000"},
              {key:'Thailand_Energy', label: 'Thailand Energy Network', linear: true, color: "#000000"},
              {key:'Vietnam_Energy', label: 'Vietnam Energy Network', linear: true, color: "#000000"}
            ]}
            tooltipLayerSources={[
              'energy_network'
            ]}
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
        <Route path="/impact">
          <Map
            map_style="impact"
            dataSources={[
              'road'
            ]}
            dataLayers={[
              {key:'road_class_1', label: 'Road Class 1', linear: true, color: "#000004"},
              {key:'road_class_2', label: 'Road Class 2', linear: true, color: "#2c115f"},
              {key:'road_class_3', label: 'Road Class 3', linear: true, color: "#721f81"},
              {key:'road_class_4', label: 'Road Class 4', linear: true, color: "#b73779"},
              {key:'road_class_5', label: 'Road Class 5', linear: true, color: "#f1605d"},
              {key:'road_class_6', label: 'Road Class 6', linear: true, color: "#feb078"}

            ]}
            tooltipLayerSources={[
              'road'
            ]}
            />
        </Route>
        <Route path="/risk">
          <Map
            map_style="risk"
            dataSources={[
              'road'
            ]}
            dataLayers={[
              {key:'road_class_1', label: 'Road Class 1', linear: true, color: "#000004"},
              {key:'road_class_2', label: 'Road Class 2', linear: true, color: "#2c115f"},
              {key:'road_class_3', label: 'Road Class 3', linear: true, color: "#721f81"},
              {key:'road_class_4', label: 'Road Class 4', linear: true, color: "#b73779"},
              {key:'road_class_5', label: 'Road Class 5', linear: true, color: "#f1605d"},
              {key:'road_class_6', label: 'Road Class 6', linear: true, color: "#feb078"}
            ]}
            tooltipLayerSources={[
              'road'
            ]}
            />
        </Route>
        <Route path="/adaptation">
          <Map
            map_style="adaptation"
            dataSources={[
              'road'
            ]}
            dataLayers={[
              {key:'road_class_1', label: 'Road Class 1', linear: true, color: "#000004"},
              {key:'road_class_2', label: 'Road Class 2', linear: true, color: "#2c115f"},
              {key:'road_class_3', label: 'Road Class 3', linear: true, color: "#721f81"},
              {key:'road_class_4', label: 'Road Class 4', linear: true, color: "#b73779"},
              {key:'road_class_5', label: 'Road Class 5', linear: true, color: "#f1605d"},
              {key:'road_class_6', label: 'Road Class 6', linear: true, color: "#feb078"}
            ]}
            tooltipLayerSources={[
              'road'
            ]}
            />
        </Route>
        <Route path="/">
          <div className="jumbotron welcome-float">
            <h1 className="h1">Southeast Asia Transport Risk Platform: Prototype Results Inquirer</h1>

            <p className="lead">
              This protoype tool presents data from the World Bank's Transport Risk Study
              of Vietnam within the risk visualisation tool developed as part of the World Bank's
              Argentia Transport Risk Study. Both studies have been undertaken for the World Bank by
              Oxford Infrastructure Analytics Ltd.
            </p>


            <p className="lead">
              The purpose of the prototype is to illustrate the type and nature of simulation results and network
              data currently available in the SEA region, and how - using the existing toolchain - this data might
              be accessed and interrogated at a National or Super-National scale.
            </p>


            <p className="lead">
              The modelling and analysis presented here aim to support decision-making by identifying
              spatial criticailities, risks, and the performance of adaptation options under current and
              future fluvial flooding outlooks. It comprises a network flow model, generation of failure
              scenarios, economic impact assessment, and cost-benefit analysis of adaptation options.
            </p>

            <p className="lead">
              The concepts and model results presented here are documented in the study report:
            </p>

            <p className="lead">
              Pant, R., Koks, E.E., Paltan, H., Russell, T., &amp; Hall, J.W. (2019). Argentina – Transport risk analysis.
              Final Report, Oxford Infrastructure Analytics Ltd., Oxford, UK. (Available by request from World Bank)
            </p>


            <p className="lead">
              The tool being used to visualize the model outputs was created and documented
              here:
            </p>
            <p className="lead">
              <a href="https://github.com/oi-analytics/argentina-transport" target="blank">Argentina Transport Study</a>
            </p>

            <p className="lead">
              <a href="https://argentina-transport-risk-analysis.readthedocs.io/en/latest/?badge=latest" target="blank">ReadTheDocs resources</a>
            </p>


            <p className="lead">
              The outputs visualized here were generated from a model created and documented
              here:
            </p>

            <p className="lead">
              <a href="https://github.com/oi-analytics/vietnam-transport" target="blank">Vietnam Transport Study</a><br></br>
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
