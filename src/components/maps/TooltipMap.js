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
        lng: lng,
        lat: lat,
        zoom: map.getZoom()
      });
    });

    map.on('mousemove', (e) => {
      const features = map.queryRenderedFeatures(e.point);
      tooltip.setLngLat(e.lngLat);

      let selectedFeatures = features.filter(features => features['source'] === 'flood')
      map.getCanvas().style.cursor = selectedFeatures.length ? 'pointer' : '';

      this.setTooltip(selectedFeatures);
    });
  }

  render() {
    const { lng, lat, zoom } = this.state;

    return (
      <div>
        <div className="custom-map-control top-left">
          <div>{`Longitude: ${lng.toFixed(2)} Latitude: ${lat.toFixed(2)} Zoom: ${zoom.toFixed(0)}`}</div>
        </div>
        <div ref={el => this.mapContainer = el} className="map" />
      </div>
    );
  }
}

TooltipMap.propTypes = {
    style: PropTypes.string.isRequired,
    tooltipLayerSources: PropTypes.array
}

export default TooltipMap
