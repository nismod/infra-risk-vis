import React from 'react'
import PropTypes from 'prop-types'

import mapboxgl from 'mapbox-gl'

class HighlightMap extends React.Component {
  tooltipContainer;

  constructor(props) {
    super(props);
    this.state = {
      lng: -56,
      lat: -30,
      zoom: 3
    };
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
      map.getCanvas().style.cursor = ''
      for (var i in features) {
        if (this.props.enabledFeatures.includes(features[i]['sourceLayer'])) {
          map.getCanvas().style.cursor = 'pointer'
        }
      }
    });

    map.on('click', (e) => {
      const features = map.queryRenderedFeatures(e.point);
      for (var i in features) {
        if (this.props.enabledFeatures.includes(features[i]['sourceLayer'])) {
          if (typeof map.getLayer('featureHighlight') !== "undefined" ){
              map.removeLayer('featureHighlight')
              map.removeSource('featureHighlight');
          }

          map.addSource('featureHighlight', {
              "type":"geojson",
              "data": features[i].toJSON()
          });
          map.addLayer({
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
          this.props.onHighlight(JSON.stringify(features[i]))
        }
      }
    })
  }

  render() {
    const { lng, lat, zoom } = this.state;

    return (
      <div>
        <div className="custom-map-control top-left">
          <div>{`Longitude: ${lng} Latitude: ${lat} Zoom: ${zoom}`}</div>
        </div>
        <div ref={el => this.mapContainer = el} className="absolute top right left bottom" />
      </div>
    );
  }
}

HighlightMap.propTypes = {
    style: PropTypes.string.isRequired,
    enabledFeatures: PropTypes.array.isRequired,
    onHighlight: PropTypes.func.isRequired
}

export default HighlightMap
