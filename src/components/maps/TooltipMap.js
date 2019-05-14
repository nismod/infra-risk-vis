import React from 'react'
import ReactDOM from 'react-dom'
import PropTypes from 'prop-types'

import mapboxgl from 'mapbox-gl'
import Tooltip from '../attributes/Tooltip'

class TooltipMap extends React.Component {
  tooltipContainer;

  constructor(props) {
    super(props);
    this.state = {
      lng: -56,
      lat: -30,
      zoom: 3
    };
  }

  setTooltip(features) {
    if (features.length) {
      ReactDOM.render(
        React.createElement(
          Tooltip, {
            features
          }
        ),
        this.tooltipContainer
      );
    } else {
      ReactDOM.render(
        React.createElement(
          'div'
        ),
        this.tooltipContainer
      );
    }
  }

  componentDidMount() {

    // Container to put React generated content in.
    this.tooltipContainer = document.createElement('div');

    const { lng, lat, zoom } = this.state;

    const map = new mapboxgl.Map({
      container: this.mapContainer,
      style: this.props.style,
      center: [lng, lat],
      zoom
    });

    const tooltip = new mapboxgl.Marker(this.tooltipContainer, {
      offset: [-120, 0]
    }).setLngLat([0,0]).addTo(map);
    
    map.on('move', () => {
      const { lng, lat } = map.getCenter();
      this.setState({
        lng: lng.toFixed(4),
        lat: lat.toFixed(4),
        zoom: map.getZoom().toFixed(2)
      });
    });

    map.on('mousemove', (e) => {
      const features = map.queryRenderedFeatures(e.point);
      tooltip.setLngLat(e.lngLat);
      
      let selectedFeatures = features.filter(features => features['source'] == 'flood')
      map.getCanvas().style.cursor = selectedFeatures.length ? 'pointer' : '';

      this.setTooltip(selectedFeatures);
    });
  }

  render() {
    const { lng, lat, zoom } = this.state;

    return (
      <div>
        <div className="absolute top left mt12 ml12 bg-darken75 color-white z1 py6 px12 round txt-s">
          <div className="mt6 mb12">
            <div className="txt-bold">Scenario</div>
            <div className="form-check">
              <input className="form-check-input" type="radio" name="scenarioRadios" id="scenarioRadios1" value="option1"/>
              <label className="form-check-label" for="scenarioRadios1">
                Baseline
              </label>
            </div>
            <div className="form-check">
              <input className="form-check-input" type="radio" name="scenarioRadios" id="scenarioRadios2" value="option2"/>
              <label className="form-check-label" for="scenarioRadios2">
                Low
              </label>
            </div>
            <div className="form-check">
              <input className="form-check-input" type="radio" name="scenarioRadios" id="scenarioRadios3" value="option3"/>
              <label className="form-check-label" for="scenarioRadios3">
                High
              </label>
            </div>

            <br/>

            <div className="txt-bold">Flood Type</div>
            <div className="form-check">
              <input className="form-check-input" type="radio" name="floodtypeRadios" id="floodtypeRadios1" value="option1"/>
              <label className="form-check-label" for="floodtypeRadios1">
                Fluvial
              </label>
            </div>
            <div className="form-check">
              <input className="form-check-input" type="radio" name="floodtypeRadios" id="floodtypeRadios2" value="option2"/>
              <label className="form-check-label" for="floodtypeRadios2">
                Pluvial
              </label>
            </div>

          </div>
        </div>

        <div className="absolute top right mt12 mr12 bg-darken75 color-white z1 py6 px12 round txt-s">
          <div className="mt6 mb12">
            <div className="txt-bold">Flood Level</div>
            <div className="form-check">
              <input className="form-check-input" type="checkbox" value="" id="floodlevel_1m-2m"/>
              <label className="form-check-label" for="floodlevel_1m-2m">
                1m-2m
              </label>
            </div>
            <div className="form-check">
              <input className="form-check-input" type="checkbox" value="" id="floodlevel_2m-3m"/>
              <label className="form-check-label" for="floodlevel_2m-3m">
                2m-3m
              </label>
            </div>
            <div className="form-check">
              <input className="form-check-input" type="checkbox" value="" id="floodlevel_3m-4m"/>
              <label className="form-check-label" for="floodlevel_3m-4m">
                3m-4m
              </label>
            </div>
            <div className="form-check">
              <input className="form-check-input" type="checkbox" value="" id="floodlevel_4m-999m"/>
              <label className="form-check-label" for="floodlevel_4m-999m">
                >4m
              </label>
            </div>
          </div>
        </div>



        <div className="absolute bottom right mb12 mr12 bg-darken75 color-white z1 py6 px12 round-full txt-s txt-bold">
          <div>{`Longitude: ${lng} Latitude: ${lat} Zoom: ${zoom}`}</div>
        </div>
        <div ref={el => this.mapContainer = el} className="absolute top right left bottom" />
      </div>
    );
  }
}

TooltipMap.propTypes = {
    style: PropTypes.string.isRequired,
    tooltipLayerSources: PropTypes.array
}  

export default TooltipMap