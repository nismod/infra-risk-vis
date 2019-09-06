import React, { Fragment } from 'react'
import ReactDOM from 'react-dom'
import PropTypes from 'prop-types'
import mapboxgl from 'mapbox-gl'

import PositionControl from './PositionControl'
import Tooltip from './Tooltip'
import FeatureSidebar from './FeatureSidebar'
import FloodHelp from './FloodHelp'
import FloodControl from './FloodControl'
import NetworkControl from './NetworkControl';

class Map extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      lng: props.lng || -56,
      lat: props.lat || -30,
      zoom: props.zoom || 4,
      selectedFeature: undefined,
      scenario: 'baseline',
      floodtype: 'fluvial',
      floodlevel: {
        _50cm1m: true,
        _1m2m: true,
        _2m3m: true,
        _3m4m: true,
        _4m999m: true
      },
      showFloodHelp: false,
      duration: 30,
      growth_rate_percentage: 2.8
    }
    this.map = undefined
    this.tooltipContainer = undefined

    this.onLayerVisChange = this.onLayerVisChange.bind(this)
    this.setScenario = this.setScenario.bind(this)
    this.setFloodType = this.setFloodType.bind(this)
    this.setFloodLevel = this.setFloodLevel.bind(this)
    this.setMap = this.setMap.bind(this)
    this.toggleFloodHelp = this.toggleFloodHelp.bind(this)
    this.updateBCR = this.updateBCR.bind(this)
  }

  setScenario(scenario) {
    this.setState({
      scenario: scenario
    })
    this.setMap(scenario, this.state.floodtype, this.state.floodlevel)
  }

  setFloodType(floodtype) {
    this.setState({
      floodtype: floodtype
    })
    this.setMap(this.state.scenario, floodtype, this.state.floodlevel)
  }

  setFloodLevel(level, value) {

    let floodlevel = Object.assign({}, this.state.floodlevel);
    floodlevel[level] = value

    this.setState({
      floodlevel: floodlevel
    })

    this.setMap(this.state.scenario, this.state.floodtype, floodlevel)
  }

  setMap(scenario, floodtype, floodlevel) {
    var flood_layers = ['50cm1m', '1m2m', '2m3m', '3m4m', '4m999m']
    var flood_layer_colors = {
      '4m999m': "#072f5f",
      '3m4m': "#1261a0",
      '2m3m': "#3895d3",
      '1m2m': "#58cced",
      '50cm1m': "#ffffff"
    }

    for (var i in flood_layers) {
      var mapLayer = this.map.getLayer('flood_' + flood_layers[i]);

      if(typeof mapLayer !== 'undefined') {
        this.map.removeLayer('flood_' + flood_layers[i]);
      }

      // insert before (under) road/rail/bridges/air/water/labels
      let before_layer_id;

      switch (this.props.map_style) {
        case 'flood':
          before_layer_id = 'country_labels';
          break;
        case 'roads':
          before_layer_id = 'road_rural';
          break;
        case 'rail':
          before_layer_id = 'rail';
          break;
        case 'airwater':
          before_layer_id = 'water';
          break;
        case 'adaptation':
          before_layer_id = 'bridges';
          break;
        case 'overview':
            before_layer_id = 'road_rural';
            break;
        default:
          before_layer_id = 'country_labels';
      }

      if (floodlevel['_' + flood_layers[i]]) {
        this.map.addLayer(
          {
            "id": "flood_" + flood_layers[i],
            "type": "fill",
            "source": "flood",
            "source-layer": scenario + '_' + floodtype + '_1in1000_' + flood_layers[i],
            "paint": {
              "fill-color": flood_layer_colors[flood_layers[i]]
            }
          },
          before_layer_id
        );
      }
    }
  }

  setTooltip(features) {
    ReactDOM.render(
      React.createElement(Tooltip, {features: features}),
      this.tooltipContainer
    );
  }

  updateBCR(data) {
    if (this.props.map_style !== 'adaptation'){
      return
    }
    const { duration, discount_growth, discount_norm, growth_rate_percentage } = data;

    const dn = discount_norm;
    const ddg = duration * discount_growth;

    const calc = [
      ">=",
      [
        "max",
        [
          "/",
          ["+",["*",["get", "baseline_ead"],dn],["*",["get", "baseline_min_eael_per_day"],ddg]],
          ["get", "baseline_tot_adap_cost"]
        ],
        [
          "/",
          ["+",["*",["get", "baseline_ead"],dn],["*",["get", "baseline_max_eael_per_day"],ddg]],
          ["get", "baseline_tot_adap_cost"]
        ],
        [
          "/",
          ["+",["*",["get", "future_med_ead"],dn],["*",["get", "future_med_min_eael_per_day"],ddg]],
          ["get", "future_med_tot_adap_cost"]
        ],
        [
          "/",
          ["+",["*",["get", "future_med_ead"],dn],["*",["get", "future_med_max_eael_per_day"],ddg]],
          ["get", "future_med_tot_adap_cost"]
        ],
        [
          "/",
          ["+",["*",["get", "future_high_ead"],dn],["*",["get", "future_high_min_eael_per_day"],ddg]],
          ["get", "future_high_tot_adap_cost"]
        ],
        [
          "/",
          ["+",["*",["get", "future_high_ead"],dn],["*",["get", "future_high_max_eael_per_day"],ddg]],
          ["get", "future_high_tot_adap_cost"]
        ]
      ],
      1
    ];

    const circle_paint_circle_color = [
      "case",
      [
        "all",
        ["has", "baseline_tot_adap_cost"],
        ["has", "future_med_tot_adap_cost"],
        ["has", "future_high_tot_adap_cost"]
      ],
      [
        "case",
        calc,
        "#860403",
        "#efefef"
      ],
      "#efefef"
    ];
    this.map.setPaintProperty('bridges', 'circle-color', circle_paint_circle_color);

    const line_paint_line_color_national = [
      "case",
      [
        "all",
        ["has", "baseline_tot_adap_cost"],
        ["has", "future_med_tot_adap_cost"],
        ["has", "future_high_tot_adap_cost"]
      ],
      [
        "case",
        calc,
        "#ba0f03",
        "#e2e2e2"
      ],
      "#e2e2e2"
    ];
    const line_paint_line_color_province = [
      "case",
      [
        "all",
        ["has", "baseline_tot_adap_cost"],
        ["has", "future_med_tot_adap_cost"],
        ["has", "future_high_tot_adap_cost"]
      ],
      [
        "case",
        calc,
        "#e0881f",
        "#e2e2e2"
      ],
      "#e2e2e2"
    ];
    const line_paint_line_color_rural = [
      "case",
      [
        "all",
        ["has", "baseline_tot_adap_cost"],
        ["has", "future_med_tot_adap_cost"],
        ["has", "future_high_tot_adap_cost"]
      ],
      [
        "case",
        calc,
        "#03ba6b",
        "#e2e2e2"
      ],
      "#e2e2e2"
    ];
    this.map.setPaintProperty('road_national', 'line-color', line_paint_line_color_national);
    this.map.setPaintProperty('road_province', 'line-color', line_paint_line_color_province);
    this.map.setPaintProperty('road_rural', 'line-color', line_paint_line_color_rural);

    this.setState({
      duration: duration,
      growth_rate_percentage: growth_rate_percentage
    })
  }

  toggleFloodHelp() {
    this.setState({
      showFloodHelp: !this.state.showFloodHelp
    })
  }

  componentDidMount() {
    const { lng, lat, zoom } = this.state

    this.map = new mapboxgl.Map({
      container: this.mapContainer,
      style: `/styles/${this.props.map_style}/style.json`,
      center: [lng, lat],
      zoom: zoom,
      minZoom: 3,
      maxZoom: 12
    })

    var nav = new mapboxgl.NavigationControl();
    this.map.addControl(nav, 'top-right');

    var scale = new mapboxgl.ScaleControl({
        maxWidth: 80,
        unit: 'metric'
    });
    this.map.addControl(scale, 'bottom-left');

    this.map.on('move', () => {
      const { lng, lat } = this.map.getCenter()
      this.setState({
        lng: lng,
        lat: lat,
        zoom: this.map.getZoom()
      })
    })

    // Container to put React generated content in.
    this.tooltipContainer = document.createElement('div')

    const tooltip = new mapboxgl.Marker(
      this.tooltipContainer, {offset: [-120, 0]}
    ).setLngLat(
      [0,0]
    ).addTo(
      this.map
    )

    this.map.on('mousemove', (e) => {
      const features = this.map.queryRenderedFeatures(e.point);

      const clickableFeatures = features.filter(
        f => this.props.dataSources.includes(f.source)
      )

      const tooltipFeatures = features.filter(
        f => this.props.tooltipLayerSources.includes(f.source)
      )

      this.map.getCanvas().style.cursor = (
        clickableFeatures.length || tooltipFeatures.length
      )? 'pointer' : '';

      tooltip.setLngLat(e.lngLat);
      this.setTooltip(tooltipFeatures);
    });

    this.map.on('click', (e) => {
      // remove current highlight
      if (typeof this.map.getLayer('featureHighlight') !== "undefined" ) {
        this.map.removeLayer('featureHighlight');
        this.map.removeSource('featureHighlight');
      }

      const features = this.map.queryRenderedFeatures(e.point);
      const clickableFeatures = features.filter(
        f => this.props.dataSources.includes(f.source)
      )

      const feature = (clickableFeatures.length)?
        clickableFeatures[0]
        : undefined;

      if (feature) {
        // add highlight layer
        this.map.addSource('featureHighlight', {
            "type":"geojson",
            "data": feature.toJSON()
        });

        if (feature.layer.type === 'line') {
          this.map.addLayer({
              "id": "featureHighlight",
              "type": "line",
              "source": "featureHighlight",
              "layout": {
                "line-join": "round",
                "line-cap": "round"
              },
              "paint": {
                "line-color": "yellow",
                "line-width": 8
              }
          });
        }
        if (feature.layer.type === 'circle') {
          this.map.addLayer({
            "id": "featureHighlight",
            "type": "circle",
            "source": "featureHighlight",
            "paint": {
              "circle-color": "yellow",
              "circle-radius": 10
            }
        });
        }
      }

      this.setState({
        selectedFeature: feature
      })
    })
  }

  onLayerVisChange(e) {
    const layer = e.target.dataset.layer
    if (e.target.checked) {
      this.map.setLayoutProperty(layer, 'visibility', 'visible')
    } else {
      this.map.setLayoutProperty(layer, 'visibility', 'none')
    }
  }

  render() {
    const { lng, lat, zoom, selectedFeature } = this.state
    const { dataLayers } = this.props

    return (
      <Fragment>
        <div className="custom-map-control top-left">
          <h3 className="h4">Select layers</h3>
          {
            (dataLayers.length)?
              <NetworkControl
                onLayerVisChange={this.onLayerVisChange}
                dataLayers={dataLayers}
              />
            : null
          }
          {
            (this.props.map_style === 'risk')?
              <small>
                Feature size indicates maximum expected annual damages plus maximum expected annual
                losses for a 30-day disruption
              </small> : null
          }
          {
            (this.props.map_style === 'impact')?
              <small>
                Feature size indicates maximum total economic impact
              </small> : null
          }
          {
            (this.props.map_style === 'roads')?
              <small>
                Feature size indicates maximum freight flows
              </small> : null
          }
          {
            (this.props.map_style === 'rail')?
              <small>
                Feature size indicates maximum freight flows
              </small> : null
          }
          {
            (this.props.map_style === 'airwater')?
              <small>
                Feature size indicates maximum freight flows (ports) or passenger flows (airports)
              </small> : null
          }
          {
            (this.props.tooltipLayerSources.includes('flood'))?
              <Fragment>
                <FloodControl
                  setScenario={this.setScenario}
                  setFloodType={this.setFloodType}
                  setFloodLevel={this.setFloodLevel}
                  />
                <a href="#flood-help" onClick={this.toggleFloodHelp}>
                { (this.state.showFloodHelp)? 'Hide info' : 'More info' }
                </a>
              </Fragment>
            : null
          }
        </div>

        <FeatureSidebar
          feature={selectedFeature}
          updateBCR={this.updateBCR}
          duration={this.state.duration}
          growth_rate_percentage={this.state.growth_rate_percentage}
          />
        { (this.state.showFloodHelp)? <FloodHelp /> : null }
        <PositionControl lat={lat} lng={lng} zoom={zoom} />
        <div ref={el => this.mapContainer = el} className="map" />
      </Fragment>
    );
  }
}

Map.propTypes = {
  map_style: PropTypes.string.isRequired,
  toggleableLayerIds: PropTypes.array,
  clickableLayerAttributes: PropTypes.object
}

export default Map
