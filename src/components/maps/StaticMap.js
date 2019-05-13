import React from 'react'
import PropTypes from 'prop-types'

import mapboxgl from 'mapbox-gl'

class StaticMap extends React.Component {

  constructor(props) {
    super(props);
    this.state = {
      lng: -56,
      lat: -30,
      zoom: 3,
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
      style: this.props.style,
      center: [lng, lat],
      zoom
    })

    this.map.on('move', () => {
      const { lng, lat } = this.map.getCenter()
      this.setState({
        lng: lng.toFixed(4),
        lat: lat.toFixed(4),
        zoom: this.map.getZoom().toFixed(2)
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
      <div>
        <div className="absolute top left mt12 ml12 bg-darken75 color-white z1 py6 px12 round txt-s">
          <div className="mt6 mb12">
            <div className="txt-bold">Layers</div>
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
        </div>

        <div className="absolute top right mt12 mr12 bg-darken75 color-white z1 py6 px12 round txt-s">
          <div className="mt6 mb12">
            <div className="txt-bold">Selected Feature</div>
            <div>{selectedFeature['_header']}</div>
            {
              Object.keys(selectedFeature).map(i => {
                if (!i.startsWith('_')) {
                  return (
                    <div>
                      {i}: {selectedFeature[i]}
                    </div>
                  )
                } else {
                  return(
                    null
                  )
                }
              })
            }

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

StaticMap.propTypes = {
  style: PropTypes.string.isRequired,
  toggleableLayerIds: PropTypes.array,
  clickableLayerAttributes: PropTypes.object
}

export default StaticMap
