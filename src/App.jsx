import React, { Fragment, useState } from 'react';
import { Route, BrowserRouter as Router, Switch } from 'react-router-dom';

import Nav from './Nav';
import Map from './Map';
import PageIntro from './PageIntro';
import RegionSummary from './RegionSummary';

import './index.css';
import 'bootstrap/dist/css/bootstrap.min.css';
import 'mapbox-gl/dist/mapbox-gl.css';

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
              {key: 'electricity', label: 'Power Grid', linear: true, color: "#eca926"},
              {key: 'rail', label: 'Railways', linear: true, color: "#444"},
              {key: 'trunk', label: 'Trunk Roads', linear: true, color: "#941339"},
              {key: 'motorway', label: 'Motorways', linear: true, color: "#941339"},
              {key: 'primary', label: 'Primary Roads', linear: true, color: "#cb3e4e"},
              {key: 'secondary', label: 'Secondary Roads', linear: true, color: "#8471a8"},
              {key: 'roads_other', label: 'Tertiary and Other Roads', linear: true, color: "#b2afaa"},

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
              {key: 'trunk', label: 'Trunk Roads', linear: true, color: "#b2afaa"},
              {key: 'motorway', label: 'Motorways', linear: true, color: "#b2afaa"},
              {key: 'primary', label: 'Primary Roads', linear: true, color: "#b2afaa"},
              {key: 'secondary', label: 'Secondary Roads', linear: true, color: "#b2afaa"},
              {key: 'roads_other', label: 'Tertiary and Other Roads', linear: true, color: "#b2afaa"},
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
              {key: 'rail', label: 'Railways', linear: true, color: "#444"},
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
              {key: 'electricity', label: 'Power Grid', linear: true, color: "#eca926"},
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
              'road',
              'flood'
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
          <PageIntro />
        </Route>
      </Switch>
      </main>
    </Fragment>
  </Router>
)
}

export default App;
