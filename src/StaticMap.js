import React, { Fragment } from 'react'
import PropTypes from 'prop-types'

import mapboxgl from 'mapbox-gl'

class StaticMap extends React.Component {

  constructor(props) {
    super(props);
    this.state = {
      lng: -56,
      lat: -30,
      zoom: 4,
      selectedFeature: {},
    }
    this.map = undefined

    this.onLayerVisClick = this.onLayerVisClick.bind(this)
    this.onFeatureHighlight = this.onFeatureHighlight.bind(this);
  }

  componentDidMount() {
    const { lng, lat, zoom } = this.state

    this.map = new mapboxgl.Map({
      container: this.mapContainer,
      style: this.props.map_style,
      center: [lng, lat],
      zoom: zoom,
      minZoom: 3,
      maxZoom: 20
    })

    var nav = new mapboxgl.NavigationControl();
    this.map.addControl(nav, 'top-right');

    this.map.on('move', () => {
      const { lng, lat } = this.map.getCenter()
      this.setState({
        lng: lng,
        lat: lat,
        zoom: this.map.getZoom()
      })
    })

    this.map.on('mousemove', (e) => {
      const features = this.map.queryRenderedFeatures(e.point);
      this.map.getCanvas().style.cursor = ''
      for (var i in features) {
        if (Object.keys(this.props.clickableLayerAttributes).includes(features[i]['layer']['id'])) {
          this.map.getCanvas().style.cursor = 'pointer'
        }
      }
    });

    this.map.on('click', (e) => {
      const features = this.map.queryRenderedFeatures(e.point);
      for (var i in features) {
        if (Object.keys(this.props.clickableLayerAttributes).includes(features[i]['layer']['id'])) {
          if (typeof this.map.getLayer('featureHighlight') !== "undefined" ){
              this.map.removeLayer('featureHighlight')
              this.map.removeSource('featureHighlight');
          }

          this.map.addSource('featureHighlight', {
              "type":"geojson",
              "data": features[i].toJSON()
          });
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
          this.onFeatureHighlight(features[i])
        }
      }
    })
  }

  onLayerVisClick(layer) {
    var visibility = this.map.getLayoutProperty(layer, 'visibility')

    if (visibility === 'visible') {
      this.map.setLayoutProperty(layer, 'visibility', 'none')
      this.className = ''
    } else {
      this.className = 'active'
      this.map.setLayoutProperty(layer, 'visibility', 'visible')
    }
  }

  onFeatureHighlight(feature) {
    const { clickableLayerAttributes } = this.props

    let attributes_names = clickableLayerAttributes[feature['layer']['id']]
    let dfeature = {}

    let attributes = Object.keys(attributes_names)
    for (var i in attributes) {

      if (attributes[i].startsWith('_')) {
        dfeature[attributes[i]] = attributes_names[attributes[i]]
      } else {
        dfeature[attributes_names[attributes[i]]] = feature['properties'][attributes[i]]
      }
    }

    this.setState({
      selectedFeature: dfeature
    })
  }

  render() {
    const { lng, lat, zoom, selectedFeature } = this.state
    const { toggleableLayerIds } = this.props

    return (
      <Fragment>
        <div className="custom-map-control top-left">
            <h4 className="h5">Layers</h4>
            {
              toggleableLayerIds.map(layer => {
                return (
                  <div className="form-check" key={'toggleLayer' + layer} >
                    <input className="form-check-input"
                      type="checkbox"
                      defaultChecked={true}
                      id={'toggleLayerCheckbox' + layer}
                      onClick={(e) => this.onLayerVisClick(layer)}/>
                    <label className="form-check-label" htmlFor={'toggleLayerCheckbox' + layer}>
                      {layer}
                    </label>
                  </div>
                )
              })
            }
        </div>

        <div className="custom-map-control top-right">
            <h4 className="h5">Selected Feature</h4>
            <div>{selectedFeature['_header']}</div>
            {
              Object.keys(selectedFeature).map(i => {
                return (i.startsWith('_'))? null : (
                    <div key={i}>
                      {i}: {selectedFeature[i]}
                    </div>
                  )
              })
            }
        </div>

        <div className="custom-map-control bottom-right">
          <div>{`Longitude: ${lng.toFixed(2)} Latitude: ${lat.toFixed(2)} Zoom: ${zoom.toFixed(0)}`}</div>
        </div>
        <div ref={el => this.mapContainer = el} className="map" />
      </Fragment>
    );
  }
}

StaticMap.propTypes = {
  map_style: PropTypes.string.isRequired,
  toggleableLayerIds: PropTypes.array,
  clickableLayerAttributes: PropTypes.object
}

export default StaticMap
