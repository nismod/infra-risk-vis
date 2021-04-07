import React, { Fragment, useState } from 'react';
import { Route, BrowserRouter as Router, Switch } from 'react-router-dom';

import Nav from './Nav'
import Map from './Map'

import './index.css';
import 'bootstrap/dist/css/bootstrap.min.css'
import 'mapbox-gl/dist/mapbox-gl.css'
import RegionSummary from './RegionSummary';

const App = () => {
  const [region, setRegion] = useState(undefined);
  return (
  <Router>
    <Fragment>
      <Route path="/" component={Nav}/>
      <main className="map-height">
      <Switch>
        <Route path="/overview">
          <Map
            map_style="overview"
            dataSources={[
              'roads',
              'rail',
              'electricity'
            ]}
            dataLayers={[
              {key: 'KHM_Electricity', label: 'Cambodia Grid', linear: true, color: "#dab540"},
              {key: 'IDN_Electricity', label: 'Indonesia Grid', linear: true, color: "#dab540"},
              {key: 'LAO_Electricity', label: 'Laos Grid', linear: true, color: "#dab540"},
              {key: 'MMR_Electricity', label: 'Myanmar Grid', linear: true, color: "#dab540"},
              {key: 'PHL_Electricity', label: 'Philippines Grid', linear: true, color: "#dab540"},
              {key: 'THA_Electricity', label: 'Thailand Grid', linear: true, color: "#dab540"},
              {key: 'VNM_Electricity', label: 'Vietnam Grid', linear: true, color: "#dab540"},
              {key: 'KHM_Rail', label: 'Cambodia Rail', linear: true, color: "#444"},
              {key: 'IDN_Rail', label: 'Indonesia Rail', linear: true, color: "#444"},
              {key: 'LAO_Rail', label: 'Laos Rail', linear: true, color: "#444"},
              {key: 'MMR_Rail', label: 'Myanmar Rail', linear: true, color: "#444"},
              {key: 'PHL_Rail', label: 'Philippines Rail', linear: true, color: "#444"},
              {key: 'THA_Rail', label: 'Thailand Rail', linear: true, color: "#444"},
              {key: 'VNM_Rail', label: 'Vietnam Rail', linear: true, color: "#444"},
              {key: 'KHM_main', label: 'Cambodia Road', linear: true, color: "#b2afaa"},
              {key: 'IDN_main', label: 'Indonesia Road', linear: true, color: "#b2afaa"},
              {key: 'LAO_main', label: 'Laos Road', linear: true, color: "#b2afaa"},
              {key: 'MMR_main', label: 'Myanmar Road', linear: true, color: "#b2afaa"},
              {key: 'PHL_main', label: 'Philippines Road', linear: true, color: "#b2afaa"},
              {key: 'THA_main', label: 'Thailand Road', linear: true, color: "#b2afaa"},
              {key: 'VNM_main', label: 'Vietnam Road', linear: true, color: "#b2afaa"}

            ]}
            tooltipLayerSources={[
              'roads',
              'rail',
              'electricity'
            ]}
            />
        </Route>
        <Route path="/roads">
          <Map
            map_style="roads"
            dataSources={[
              'roads'
            ]}
            dataLayers={[
              {key: 'KHM_main', label: 'Cambodia Major Roads', linear: true, color: "#e48d14"},
              {key: 'IDN_main', label: 'Indonesia Major Roads', linear: true, color: "#e48d14"},
              {key: 'LAO_main', label: 'Laos Major Roads', linear: true, color: "#e48d14"},
              {key: 'MMR_main', label: 'Myanmar Major Roads', linear: true, color: "#e48d14"},
              {key: 'PHL_main', label: 'Philippines Major Roads', linear: true, color: "#e48d14"},
              {key: 'THA_main', label: 'Thailand Major Roads', linear: true, color: "#e48d14"},
              {key: 'VNM_main', label: 'Vietnam Major Roads', linear: true, color: "#e48d14"},
              {key: 'KHM_other', label: 'Cambodia Minor Roads', linear: true, color: "#e48d14"},
              {key: 'IDN_other', label: 'Indonesia Minor Roads', linear: true, color: "#b2afaa"},
              {key: 'LAO_other', label: 'Laos Minor Roads', linear: true, color: "#b2afaa"},
              {key: 'MMR_other', label: 'Myanmar Minor Roads', linear: true, color: "#b2afaa"},
              {key: 'PHL_other', label: 'Philippines Minor Roads', linear: true, color: "#b2afaa"},
              {key: 'THA_other', label: 'Thailand Minor Roads', linear: true, color: "#b2afaa"},
              {key: 'VNM_other', label: 'Vietnam Minor Roads', linear: true, color: "#b2afaa"}

            ]}
            tooltipLayerSources={[
              'roads'
            ]}
            />
        </Route>
        <Route path="/rail">
          <Map
            map_style="rail"
            dataSources={[
              'rail'
            ]}
            dataLayers={[
              {key: 'KHM_Rail', label: 'Cambodia Rail', linear: true, color: "#444"},
              {key: 'IDN_Rail', label: 'Indonesia Rail', linear: true, color: "#444"},
              {key: 'LAO_Rail', label: 'Laos Rail', linear: true, color: "#444"},
              {key: 'MMR_Rail', label: 'Myanmar Rail', linear: true, color: "#444"},
              {key: 'PHL_Rail', label: 'Philippines Rail', linear: true, color: "#444"},
              {key: 'THA_Rail', label: 'Thailand Rail', linear: true, color: "#444"},
              {key: 'VNM_Rail', label: 'Vietnam Rail', linear: true, color: "#444"}
            ]}
            tooltipLayerSources={[
              'rail'
            ]}
            />
        </Route>
        <Route path="/energy_network">
          <Map
            map_style="electricity"
            dataSources={[
              'electricity'
            ]}
            dataLayers={[
              {key: 'KHM_Electricity', label: 'Cambodia Grid', linear: true, color: "#dab540"},
              {key: 'IDN_Electricity', label: 'Indonesia Grid', linear: true, color: "#dab540"},
              {key: 'LAO_Electricity', label: 'Laos Grid', linear: true, color: "#dab540"},
              {key: 'MMR_Electricity', label: 'Myanmar Grid', linear: true, color: "#dab540"},
              {key: 'PHL_Electricity', label: 'Philippines Grid', linear: true, color: "#dab540"},
              {key: 'THA_Electricity', label: 'Thailand Grid', linear: true, color: "#dab540"},
              {key: 'VNM_Electricity', label: 'Vietnam Grid', linear: true, color: "#dab540"},
            ]}
            tooltipLayerSources={[
              'electricity'
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
        <Route path="/summary">
          <div className="page-col-right">
            <RegionSummary region={region} />
          </div>
          <div className="page-col-left">
            <Map
              map_style="regions"
              dataSources={[
                'boundaries'
              ]}
              dataLayers={[]}
              tooltipLayerSources={[]}
              onRegionSelect={setRegion}
              />
          </div>
        </Route>
        <Route path="/">
          <article>
            <h1 className="h1">Southeast Asia Infrastructure Risk: Prototype</h1>

            <p>
              This protoype tool presents data from the World Bank's Transport Risk Study
              of Vietnam within the risk visualisation tool developed as part of the World Bank's
              Argentia Transport Risk Study. Both studies have been undertaken for the World Bank by
              Oxford Infrastructure Analytics Ltd.
            </p>


            <p>
              The purpose of the prototype is to illustrate the type and nature of simulation results and network
              data currently available in the SEA region, and how - using the existing toolchain - this data might
              be accessed and interrogated at a National or Super-National scale.
            </p>


            <p>
              The modelling and analysis presented here aim to support decision-making by identifying
              spatial criticailities, risks, and the performance of adaptation options under current and
              future fluvial flooding outlooks. It comprises a network flow model, generation of failure
              scenarios, economic impact assessment, and cost-benefit analysis of adaptation options.
            </p>

            <p>
              The concepts and model results presented here are documented in the study report:
            </p>

            <p>
              Pant, R., Koks, E.E., Paltan, H., Russell, T., &amp; Hall, J.W. (2019). Argentina – Transport risk analysis.
              Final Report, Oxford Infrastructure Analytics Ltd., Oxford, UK. (Available by request from World Bank)
            </p>


            <p>
              The tool being used to visualize the model outputs was created and documented
              here:
            </p>
            <p>
              <a href="https://github.com/oi-analytics/argentina-transport" target="blank">Argentina Transport Study</a>
            </p>

            <p>
              <a href="https://argentina-transport-risk-analysis.readthedocs.io/en/latest/?badge=latest" target="blank">ReadTheDocs resources</a>
            </p>


            <p>
              The outputs visualized here were generated from a model created and documented
              here:
            </p>

            <p>
              <a href="https://github.com/oi-analytics/vietnam-transport" target="blank">Vietnam Transport Study</a><br></br>
            </p>


            <h1 className="h1">Funding support</h1>
            <p>

            This results inquirer tool has been developed for the Government of Argentina with
            funding support from the World Bank Group and Global Facility for Disaster
            Reduction and Recovery (GFDRR).

            </p>
          </article>
        </Route>
      </Switch>
      </main>
    </Fragment>
  </Router>
)
}

export default App;
