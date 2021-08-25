import React, { createRef, Fragment, useCallback, useEffect, useRef, useState } from 'react';
import PropTypes from 'prop-types';
import mapboxgl from 'mapbox-gl';
import Drawer from '@material-ui/core/Drawer';
import Toolbar from '@material-ui/core/Toolbar';

import BackgroundControl from './BackgroundControl';
import PositionControl from './PositionControl';
import Tooltip from './Tooltip';
import FeatureSidebar from './FeatureSidebar';
import Help from './Help';
import FloodControl from './FloodControl';
import NetworkControl from './NetworkControl';
import RiskControl from './RiskControl';
import { createPortal } from 'react-dom';

const MAPBOX_KEY = 'pk.eyJ1IjoidG9tcnVzc2VsbCIsImEiOiJjaXZpMTFpdGkwMDQ1MnptcTh4ZzRzeXNsIn0.ZSvSOHSsWBQ44QNhA71M6Q';
mapboxgl.accessToken = MAPBOX_KEY;

function networkBaseLayerID(mapStyle) {
  // insert before (under) road/rail/bridges/air/water/labels
  let before_layer_id;

  switch (mapStyle) {
    case 'flood':
      before_layer_id = 'country_labels';
      break;
    case 'risk':
      before_layer_id = 'road_class_6';
      break;
    default:
      before_layer_id = 'country_labels';
  }
  return before_layer_id;
}

function drawFeature(feature, map) {
  // remove current highlight
  if (map.getLayer('featureHighlight')) {
    map.removeLayer('featureHighlight');
  }

  // add highlight layer
  const source = map.getSource('featureHighlight');
  if (source) {
    source.setData(feature.toJSON());
  } else {
    map.addSource('featureHighlight', {
      type: 'geojson',
      data: feature.toJSON(),
    });
  }

  if (feature.layer.type === 'line') {
    map.addLayer(
      {
        id: 'featureHighlight',
        type: 'line',
        source: 'featureHighlight',
        layout: {
          'line-join': 'round',
          'line-cap': 'round',
        },
        paint: {
          'line-color': 'yellow',
          'line-width': {
            base: 1,
            stops: [
              [3, 1],
              [10, 8],
              [17, 16],
            ],
          },
        },
      },
      feature.layer.id,
    );
  }
  if (feature.layer.type === 'circle') {
    map.addLayer(
      {
        id: 'featureHighlight',
        type: 'circle',
        source: 'featureHighlight',
        paint: {
          'circle-color': 'yellow',
          'circle-radius': {
            base: 1,
            stops: [
              [3, 4],
              [10, 12],
              [17, 20],
            ],
          },
        },
      },
      feature.layer.id,
    );
  }
}

const Map = ({
  lng: lngProp,
  lat: latProp,
  zoom: zoomProp,
  map_style,
  dataLayers,
  tooltipLayerSources,
  dataSources,
  onRegionSelect,
}) => {
  const [lng, setLng] = useState(lngProp ?? -77.28);
  const [lat, setLat] = useState(latProp ?? 18.14);
  const [zoom, setZoom] = useState(zoomProp ?? 8);

  const [selectedFeature, setSelectedFeature] = useState();
  const [layerVisibility, setLayerVisibility] = useState(Object.fromEntries(dataLayers.map((dl) => [dl.key, true])));

  const [background, setBackground] = useState('light');

  const [showHelp, setShowHelp] = useState(false);
  const [helpTopic, setHelpTopic] = useState();

  const [riskMetric, setRiskMetric] = useState('total');
  const [duration, setDuration] = useState(30);
  const [growth_rate_percentage, setGrowth_rate_percentage] = useState(2.8);

  const [scenario, setScenario] = useState('baseline');
  const [floodlevel, setFloodLevel] = useState({
    _1m2m: false,
    _2m3m: false,
    _3m4m: false,
    _4m999m: false,
  });
  const [floodtype, setFloodType] = useState('fluvial');

  const updateFloodLevel = useCallback(
    (level, value) => {
      setFloodLevel({ ...floodlevel, [level]: value });
    },
    [floodlevel],
  );

  const [map, setMap] = useState();
  const mapContainerRef = createRef();
  const tooltip = useRef();
  const tooltipContainerRef = useRef(document.createElement('div'));
  const [tooltipFeatures, setTooltipFeatures] = useState();

  useEffect(() => {
    const newMap = new mapboxgl.Map({
      container: mapContainerRef.current,
      style: `/styles/${map_style}/style.json`,
      center: [lng, lat],
      zoom: zoom,
      minZoom: 3,
      maxZoom: 16,
    });

    newMap.on('load', () => {
      var nav = new mapboxgl.NavigationControl();
      newMap.addControl(nav, 'top-right');

      var scale = new mapboxgl.ScaleControl({
        maxWidth: 80,
        unit: 'metric',
      });
      newMap.addControl(scale, 'bottom-left');

      // offset to match width set in tooltip-body css
      tooltip.current = new mapboxgl.Marker(tooltipContainerRef.current, { offset: [-150, 0] })
        .setLngLat([0, 0])
        .addTo(newMap);

      setMap(newMap);
    });
  }, []);

  const onMapMove = useCallback(() => {
    const { lng, lat } = map.getCenter();
    setLat(lat);
    setLng(lng);
    setZoom(map.getZoom());
  }, [map]);

  const onMapMouseMove = useCallback(
    (e) => {
      const features = map.queryRenderedFeatures(e.point);

      const clickableFeatures = features.filter((f) => dataSources.includes(f.source));

      const tooltipFeatures = features.filter((f) => tooltipLayerSources.includes(f.source));

      map.getCanvas().style.cursor = clickableFeatures.length || tooltipFeatures.length ? 'pointer' : '';

      tooltip.current.setLngLat(e.lngLat);
      setTooltipFeatures(tooltipFeatures);
    },
    [dataSources, map, tooltipLayerSources],
  );

  const onMapClick = useCallback(
    (e) => {
      const features = map.queryRenderedFeatures(e.point);
      const clickableFeatures = features.filter((f) => dataSources.includes(f.source));

      const feature = clickableFeatures?.[0];

      if (map_style === 'regions') {
        // pass region code up to App for RegionSummary to use
        onRegionSelect(feature?.properties);
      } else {
        if (feature) {
          drawFeature(feature, map);
        } else {
          // remove current highlight
          if (map.getLayer('featureHighlight')) {
            map.removeLayer('featureHighlight');
          }
        }
        setSelectedFeature(feature);
      }
    },
    [dataSources, map, map_style, onRegionSelect],
  );

  useEffect(() => {
    const handler = onMapMove;
    map?.on('move', handler);

    return () => {
      map?.off('move', handler);
    };
  }, [map, onMapMove]);

  useEffect(() => {
    const handler = onMapMouseMove;
    map?.on('mousemove', handler);

    return () => {
      map?.off('mousemove', handler);
    };
  }, [map, onMapMouseMove]);

  useEffect(() => {
    const handler = onMapClick;
    map?.on('click', handler);

    return () => {
      map?.off('click', handler);
    };
  }, [map, onMapClick]);

  useEffect(() => {
    if (!map) return;

    const flood_layers = ['1m2m', '2m3m', '3m4m', '4m999m'];
    const flood_layer_colors = {
      '4m999m': '#072f5f',
      '3m4m': '#1261a0',
      '2m3m': '#3895d3',
      '1m2m': '#58cced',
    };

    for (const layer of flood_layers) {
      const layerKey = 'flood_' + layer;
      const mapLayer = map.getLayer(layerKey);

      if (typeof mapLayer !== 'undefined') {
        map.removeLayer(layerKey);
      }

      if (floodlevel['_' + layer]) {
        map.addLayer(
          {
            id: layerKey,
            type: 'fill',
            source: 'flood',
            'source-layer': `${scenario}_${floodtype}_1in1000_${layer}`,
            paint: {
              'fill-color': flood_layer_colors[layer],
            },
          },
          networkBaseLayerID(map_style),
        );
      }
    }
  }, [scenario, floodlevel, floodtype, map_style, map]);

  useEffect(() => {
    let calc;

    if (map_style !== 'electricity') {
      // EAD has EAEL baked in to the numbers in road and rail
      if (riskMetric === 'total') {
        // EAD subtract 335/365 of EAEL (annual) to get total
        // assuming 30-day disruption
        calc = ['-', ['get', 'EAD_max'], ['*', 0.917808219178, ['coalesce', ['get', 'EAEL'], 0]]];
      }
      if (riskMetric === 'EAD') {
        // EAD subtract EAEL to get direct only
        calc = ['-', ['get', 'EAD_max'], ['coalesce', ['get', 'EAEL'], 0]];
      }
      if (riskMetric === 'EAEL') {
        calc = ['*', 0.0821917808219, ['coalesce', ['get', 'EAEL'], 0]];
      }
    } else {
      // EAD does not include EAEL in electricity
      // EAEL here is in annual dollars, so multiply by
      // 1e-6 * (30/365)
      // to get 30-day disruption in US$m
      if (riskMetric === 'total') {
        calc = ['+', ['get', 'EAD_max'], ['*', 0.0000000821917808219, ['coalesce', ['get', 'EAEL'], 0]]];
      }
      if (riskMetric === 'EAD') {
        calc = ['get', 'EAD_max'];
      }
      if (riskMetric === 'EAEL') {
        calc = ['*', 0.0000000821917808219, ['coalesce', ['get', 'EAEL'], 0]];
      }
    }

    const paint_color = [
      'interpolate-lab',
      ['linear'],
      calc,
      0,
      ['to-color', '#b2afaa'],
      0.000000001,
      ['to-color', '#fff'],
      0.001,
      ['to-color', '#fcfcb8'],
      0.01,
      ['to-color', '#ff9c66'],
      0.1,
      ['to-color', '#d03f6f'],
      1,
      ['to-color', '#792283'],
      10,
      ['to-color', '#3f0a72'],
      100,
      ['to-color', '#151030'],
    ];

    if (map_style === 'roads') {
      map.setPaintProperty('trunk', 'line-color', paint_color);
      map.setPaintProperty('primary', 'line-color', paint_color);
      map.setPaintProperty('secondary', 'line-color', paint_color);
      map.setPaintProperty('roads_other', 'line-color', paint_color);
      map.setPaintProperty('motorway', 'line-color', paint_color);
    }
    if (map_style === 'rail' || map_style === 'electricity') {
      map.setPaintProperty(map_style, 'line-color', paint_color);
    }
  }, [map, map_style, riskMetric]);

  useEffect(() => {
    if (!map || !map_style) return;

    fetch(`/styles/${map_style}/style.json`)
      .then((response) => response.json())
      .then((data) => {
        map.setStyle(data);
        if (selectedFeature && dataSources.includes(selectedFeature.source)) {
          drawFeature(selectedFeature, map);
        }

        // TODO - understand why this is called here?
        // if (map_style === 'roads' || map_style === 'rail' || map_style === 'electricity') {
        //   this.setRiskMetric(this.state.riskMetric);
        // }
      });
  }, [dataSources, map, map_style, selectedFeature]);

  const toggleHelp = useCallback(
    (e) => {
      const newHelpTopic = e.target.dataset.helpTopic;
      const newShowHelp = !showHelp || helpTopic !== newHelpTopic;
      setShowHelp(newShowHelp);
      setHelpTopic(newHelpTopic);
    },
    [showHelp, helpTopic],
  );

  const updateBCR = useCallback(({ duration, growth_rate_percentage }) => {
    setDuration(duration);
    setGrowth_rate_percentage(growth_rate_percentage);
  }, []);

  const onBackgroundChange = useCallback(
    (e) => {
      const newBackground = e.target.value;

      map.setLayoutProperty(newBackground, 'visibility', 'visible');
      map.setLayoutProperty(background, 'visibility', 'none');
      setBackground(newBackground);
    },
    [background, map],
  );

  const onLayerVisChange = useCallback(
    (e) => {
      const layer = e.target.value;
      const value = e.target.checked;
      if (value) {
        map.setLayoutProperty(layer, 'visibility', 'visible');
      } else {
        map.setLayoutProperty(layer, 'visibility', 'none');
      }
      setLayerVisibility({ ...layerVisibility, [layer]: value });
    },
    [layerVisibility, map],
  );

  return (
    <>
      <Drawer variant="permanent">
        <Toolbar />
        <div className="drawer-contents">
          {dataLayers.length ? (
            <>
              <h2 className="h4">Select layers</h2>
              <NetworkControl
                onLayerVisChange={onLayerVisChange}
                dataLayers={dataLayers}
                layerVisibility={layerVisibility}
              />
              <BackgroundControl onBackgroundChange={onBackgroundChange} background={background} />
            </>
          ) : null}
          {map_style === 'regions' ? (
            <>
              <small>Max Total Expected Risk (EAD + EAEL for 30 day disruption, million US$)</small>
              <svg width="270" height="25" version="1.1" xmlns="http://www.w3.org/2000/svg">
                <defs>
                  <linearGradient id="summary_gradient" x1="0" x2="1" y1="0" y2="0">
                    <stop offset="0%" stopColor="#ffffff" />
                    <stop offset="12.5%" stopColor="#fee0d2" />
                    <stop offset="25%" stopColor="#fdc1a9" />
                    <stop offset="37.5%" stopColor="#fc9d7f" />
                    <stop offset="50%" stopColor="#fb7859" />
                    <stop offset="62.5%" stopColor="#f4513b" />
                    <stop offset="75%" stopColor="#de2c26" />
                    <stop offset="87.5%" stopColor="#bf161b" />
                    <stop offset="100%" stopColor="#950b13" />
                  </linearGradient>
                </defs>
                <g fill="none" fontSize="10" fontFamily="sans-serif"></g>
                <rect x="2" y="0" width="258" height="10" fill="url(#summary_gradient)" />
                <g fill="none" fontSize="10" transform="translate(2,10)" fontFamily="sans-serif" textAnchor="middle">
                  <g transform="translate(0.5,0)">
                    <line stroke="currentColor" y2="3"></line>
                    <text fill="currentColor" y="6" dy="0.71em">
                      0
                    </text>
                  </g>
                  <g transform="translate(32.5,0)">
                    <line stroke="currentColor" y2="3"></line>
                    <text fill="currentColor" y="6" dy="0.71em">
                      0.1
                    </text>
                  </g>
                  <g transform="translate(65,0)">
                    <line stroke="currentColor" y2="3"></line>
                    <text fill="currentColor" y="6" dy="0.71em">
                      0.5
                    </text>
                  </g>
                  <g transform="translate(97.5,0)">
                    <line stroke="currentColor" y2="3"></line>
                    <text fill="currentColor" y="6" dy="0.71em">
                      2.5
                    </text>
                  </g>
                  <g transform="translate(130,0)">
                    <line stroke="currentColor" y2="3"></line>
                    <text fill="currentColor" y="6" dy="0.71em">
                      10
                    </text>
                  </g>
                  <g transform="translate(162.5,0)">
                    <line stroke="currentColor" y2="3"></line>
                    <text fill="currentColor" y="6" dy="0.71em">
                      50
                    </text>
                  </g>
                  <g transform="translate(195,0)">
                    <line stroke="currentColor" y2="3"></line>
                    <text fill="currentColor" y="6" dy="0.71em">
                      250
                    </text>
                  </g>
                  <g transform="translate(227.5,0)">
                    <line stroke="currentColor" y2="3"></line>
                    <text fill="currentColor" y="6" dy="0.71em">
                      1k
                    </text>
                  </g>
                  <g transform="translate(257.5,0)">
                    <line stroke="currentColor" y2="3"></line>
                    <text fill="currentColor" y="6" dy="0.71em">
                      5k
                    </text>
                  </g>
                </g>
              </svg>
            </>
          ) : null}
          {map_style === 'risk' ? (
            <div>
              <small>
                Feature size indicates maximum expected annual damages plus maximum expected annual losses for a 30-day
                disruption
              </small>
              <span className="dot line" style={{ height: '2px', width: '24px' }}></span>&lt;1 million USD
              <br />
              <span className="dot line" style={{ height: '4px', width: '24px' }}></span>1-5 million USD
              <br />
              <span className="dot line" style={{ height: '6px', width: '24px' }}></span>5-10 million USD
              <br />
              <span className="dot line" style={{ height: '8px', width: '24px' }}></span>&gt;10 million USD
              <br />
              <a href="#help" data-help-topic="vietnam" onClick={toggleHelp}>
                {showHelp && helpTopic === 'vietnam' ? 'Hide info' : 'More info'}
              </a>
            </div>
          ) : null}
          {map_style === 'roads' || map_style === 'rail' || map_style === 'electricity' ? (
            <RiskControl setRiskMetric={setRiskMetric} riskMetric={riskMetric} />
          ) : null}
          {map_style === 'roads' ? <small>Road network data extracted from OpenStreetMap</small> : null}
          {map_style === 'rail' ? (
            <div>
              <small>Rail network data extracted from OpenStreetMap</small>
            </div>
          ) : null}
          {map_style === 'electricity' ? <small>Energy network data extracted from Gridfinder</small> : null}
          {map_style === 'hazards' ? (
            <>
              <small>Coastal flood depth (m)</small>
              <svg width="270" height="25" version="1.1" xmlns="http://www.w3.org/2000/svg">
                <defs>
                  <linearGradient id="coastal_gradient" x1="0" x2="1" y1="0" y2="0">
                    <stop offset="0%" stopColor="#9df4b0" />
                    <stop offset="100%" stopColor="#0b601d" />
                  </linearGradient>
                </defs>
                <g fill="none" fontSize="10" fontFamily="sans-serif"></g>
                <rect x="2" y="0" width="258" height="10" fill="url(#coastal_gradient)" />
                <g fill="none" fontSize="10" transform="translate(2,10)" fontFamily="sans-serif" textAnchor="middle">
                  <g transform="translate(0.5,0)">
                    <line stroke="currentColor" y2="3"></line>
                    <text fill="currentColor" y="6" dy="0.71em">
                      0
                    </text>
                  </g>
                  <g transform="translate(130,0)">
                    <line stroke="currentColor" y2="3"></line>
                    <text fill="currentColor" y="6" dy="0.71em">
                      2.5
                    </text>
                  </g>
                  <g transform="translate(257.5,0)">
                    <line stroke="currentColor" y2="3"></line>
                    <text fill="currentColor" y="6" dy="0.71em">
                      5
                    </text>
                  </g>
                </g>
              </svg>
              <small>Fluvial flood depth (m)</small>
              <svg width="270" height="25" version="1.1" xmlns="http://www.w3.org/2000/svg">
                <defs>
                  <linearGradient id="fluvial_gradient" x1="0" x2="1" y1="0" y2="0">
                    <stop offset="0%" stopColor="#58cced" />
                    <stop offset="100%" stopColor="#072f5f" />
                  </linearGradient>
                </defs>
                <g fill="none" fontSize="10" fontFamily="sans-serif"></g>
                <rect x="2" y="0" width="258" height="10" fill="url(#fluvial_gradient)" />
                <g fill="none" fontSize="10" transform="translate(2,10)" fontFamily="sans-serif" textAnchor="middle">
                  <g transform="translate(0.5,0)">
                    <line stroke="currentColor" y2="3"></line>
                    <text fill="currentColor" y="6" dy="0.71em">
                      0
                    </text>
                  </g>
                  <g transform="translate(130,0)">
                    <line stroke="currentColor" y2="3"></line>
                    <text fill="currentColor" y="6" dy="0.71em">
                      2.5
                    </text>
                  </g>
                  <g transform="translate(257.5,0)">
                    <line stroke="currentColor" y2="3"></line>
                    <text fill="currentColor" y="6" dy="0.71em">
                      5
                    </text>
                  </g>
                </g>
              </svg>
              <small>Cyclone gust speed (m/s)</small>
              <svg width="270" height="25" version="1.1" xmlns="http://www.w3.org/2000/svg">
                <defs>
                  <linearGradient id="cyclone_gradient" x1="0" x2="1" y1="0" y2="0">
                    <stop offset="0%" stopColor="#ffffff" />
                    <stop offset="50%" stopColor="#f9d5cb" />
                    <stop offset="100%" stopColor="#d44118" />
                  </linearGradient>
                </defs>
                <g fill="none" fontSize="10" fontFamily="sans-serif"></g>
                <rect x="2" y="0" width="258" height="10" fill="url(#cyclone_gradient)" />
                <g fill="none" fontSize="10" transform="translate(2,10)" fontFamily="sans-serif" textAnchor="middle">
                  <g transform="translate(0.5,0)">
                    <line stroke="currentColor" y2="3"></line>
                    <text fill="currentColor" y="6" dy="0.71em">
                      0
                    </text>
                  </g>
                  <g transform="translate(130,0)">
                    <line stroke="currentColor" y2="3"></line>
                    <text fill="currentColor" y="6" dy="0.71em">
                      25
                    </text>
                  </g>
                  <g transform="translate(257.5,0)">
                    <line stroke="currentColor" y2="3"></line>
                    <text fill="currentColor" y="6" dy="0.71em">
                      50
                    </text>
                  </g>
                </g>
              </svg>
              <a href="#help" data-help-topic="hazards" onClick={toggleHelp}>
                {showHelp && helpTopic === 'hazards' ? 'Hide info' : 'More info'}
              </a>
            </>
          ) : null}
          {tooltipLayerSources.includes('flood') ? (
            <>
              <FloodControl setScenario={setScenario} setFloodType={setFloodType} setFloodLevel={updateFloodLevel} />
              <a href="#help" data-help-topic="flood" onClick={toggleHelp}>
                {showHelp && helpTopic === 'flood' ? 'Hide info' : 'More info'}
              </a>
            </>
          ) : null}
        </div>
      </Drawer>
      <div className="map-height">
        {showHelp ? (
          <Help topic={helpTopic} />
        ) : (
          <FeatureSidebar
            feature={selectedFeature}
            updateBCR={updateBCR}
            duration={duration}
            growth_rate_percentage={growth_rate_percentage}
          />
        )}
        <PositionControl lat={lat} lng={lng} zoom={zoom} />
        <div ref={mapContainerRef} className="map" />
      </div>
      {tooltipFeatures &&
        createPortal(<Tooltip features={tooltipFeatures} map_style={map_style}></Tooltip>, tooltipContainerRef.current)}
    </>
  );
};

Map.propTypes = {
  map_style: PropTypes.string.isRequired,
  lat: PropTypes.number,
  lng: PropTypes.number,
  zoom: PropTypes.number,
  dataSources: PropTypes.array.isRequired,
  dataLayers: PropTypes.array.isRequired,
  tooltipLayerSources: PropTypes.array.isRequired,
  onRegionSelect: PropTypes.func,
};

export default Map;
