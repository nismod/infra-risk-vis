import React, { Fragment } from 'react'
import ReactDOM from 'react-dom'
import PropTypes from 'prop-types'

import mapboxgl from 'mapbox-gl'
import Tooltip from './Tooltip'
import PositionControl from './PositionControl'

class TooltipMap extends React.Component {
  tooltipContainer;

  constructor(props) {
    super(props);
    this.state = {
      lng: -56,
      lat: -30,
      zoom: 3,
      scenario: 'baseline',
      floodtype: 'fluvial',
      floodlevel: {
        _1m2m: true,
        _2m3m: true,
        _3m4m: true,
        _4m999m: true
      }
    };

    this.map = undefined
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

    var flood_layers = ['1m2m', '2m3m', '3m4m', '4m999m']
    var flood_layer_colors = {
      '1m2m': "#072f5f",
      '2m3m': "#1261a0",
      '3m4m': "#3895d3",
      '4m999m': "#58cced"
    }

    for (var i in flood_layers) {

      var mapLayer = this.map.getLayer('flood_' + flood_layers[i]);

      if(typeof mapLayer !== 'undefined') {
        this.map.removeLayer('flood_' + flood_layers[i]);
      }

      if (floodlevel['_' + flood_layers[i]]) {
        this.map.addLayer({
          "id": "flood_" + flood_layers[i],
          "type": "fill",
          "source": "flood",
          "source-layer": scenario + '_' + floodtype + '_1in500_' + flood_layers[i],
          "paint": {
            "fill-color": flood_layer_colors[flood_layers[i]]
          }
        },);
      }
    }
  }

  setTooltip(features) {
    ReactDOM.render(
      React.createElement(Tooltip, {features: features}),
      this.tooltipContainer
    );
  }

  componentDidMount() {

    // Container to put React generated content in.
    this.tooltipContainer = document.createElement('div');

    const { lng, lat, zoom } = this.state;

    this.map = new mapboxgl.Map({
      container: this.mapContainer,
      style: this.props.map_style,
      center: [lng, lat],
      zoom: zoom,
      minZoom: 3,
      maxZoom: 12
    });

    const tooltip = new mapboxgl.Marker(this.tooltipContainer, {
      offset: [-120, 0]
    }).setLngLat([0,0]).addTo(this.map);

    this.map.on('move', () => {
      const { lng, lat } = this.map.getCenter();
      this.setState({
        lng: lng,
        lat: lat,
        zoom: this.map.getZoom()
      });
    });

    this.map.on('mousemove', (e) => {
      const features = this.map.queryRenderedFeatures(e.point);
      tooltip.setLngLat(e.lngLat);

      let selectedFeatures = features.filter(features => features['source'] === 'flood')
      this.map.getCanvas().style.cursor = selectedFeatures.length ? 'pointer' : '';

      this.setTooltip(selectedFeatures);
    });
  }

  render() {
    const { lng, lat, zoom } = this.state;

    return (
      <Fragment>
        <div className="custom-map-control top-left">
            <h4 className="h5">Scenario</h4>
            <div className="form-check">
              <input className="form-check-input" defaultChecked={true} type="radio" name="scenarioRadio" value="baseline" onClick={(e) => this.setScenario(e.target.value)}/>
              <label className="form-check-label">
                Baseline
              </label>
            </div>
            <div className="form-check">
              <input className="form-check-input" type="radio" name="scenarioRadio" value="low" onClick={(e) => this.setScenario(e.target.value)}/>
              <label className="form-check-label">
                Low
              </label>
            </div>
            <div className="form-check">
              <input className="form-check-input" type="radio" name="scenarioRadio" value="high" onClick={(e) => this.setScenario(e.target.value)}/>
              <label className="form-check-label">
                High
              </label>
            </div>

            <br/>

            <h4 className="h5">Flood Type</h4>
            <div className="form-check">
              <input className="form-check-input" defaultChecked={true} type="radio" name="floodtypeRadios" value="fluvial" onClick={(e) => this.setFloodType(e.target.value)}/>
              <label className="form-check-label">
                Fluvial
              </label>
            </div>
            <div className="form-check">
              <input className="form-check-input" type="radio" name="floodtypeRadios" value="pluvial" onClick={(e) => this.setFloodType(e.target.value)}/>
              <label className="form-check-label">
                Pluvial
              </label>
            </div>
        </div>

        <div className="custom-map-control top-right">
          <h4 className="h5">Flood Level</h4>
          <div className="form-check">
            <input className="form-check-input" defaultChecked={true} type="checkbox" value="_1m2m" onClick={(e) => this.setFloodLevel(e.target.value, e.target.checked)}/>
            <label className="form-check-label">
              1m-2m
            </label>
          </div>
          <div className="form-check">
            <input className="form-check-input" defaultChecked={true} type="checkbox" value="_2m3m" onClick={(e) => this.setFloodLevel(e.target.value, e.target.checked)}/>
            <label className="form-check-label">
              2m-3m
            </label>
          </div>
          <div className="form-check">
            <input className="form-check-input" defaultChecked={true} type="checkbox" value="_3m4m" onClick={(e) => this.setFloodLevel(e.target.value, e.target.checked)}/>
            <label className="form-check-label">
              3m-4m
            </label>
          </div>
          <div className="form-check">
            <input className="form-check-input" defaultChecked={true} type="checkbox" value="_4m999m" onClick={(e) => this.setFloodLevel(e.target.value, e.target.checked)}/>
            <label className="form-check-label">
              >4m
            </label>
          </div>
        </div>

        <PositionControl lat={lat} lng={lng} zoom={zoom} />
        <div ref={el => this.mapContainer = el} className="map" />
      </Fragment>
    );
  }
}

TooltipMap.propTypes = {
    map_style: PropTypes.string.isRequired,
    tooltipLayerSources: PropTypes.array
}

export default TooltipMap
