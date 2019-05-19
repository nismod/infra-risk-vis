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
      selectedFeature: undefined
    }
    this.map = undefined

    this.onLayerVisChange = this.onLayerVisChange.bind(this)
  }

  componentDidMount() {
    const { lng, lat, zoom } = this.state

    this.map = new mapboxgl.Map({
      container: this.mapContainer,
      style: this.props.map_style,
      center: [lng, lat],
      zoom: zoom,
      minZoom: 3,
      maxZoom: 12
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
        if (this.props.dataSources.includes(features[i].source)) {
          this.map.getCanvas().style.cursor = 'pointer'
          break
        }
      }
    });

    this.map.on('click', (e) => {
      // remove current highlight
      if (typeof this.map.getLayer('featureHighlight') !== "undefined" ) {
        this.map.removeLayer('featureHighlight');
        this.map.removeSource('featureHighlight');
      }

      let feature;
      const features = this.map.queryRenderedFeatures(e.point);
      for (var i in features) {
        if (this.props.dataSources.includes(features[i].source)) {
          feature = features[i]
          break
        }
      }

      if (feature) {
        // add highlight layer
        this.map.addSource('featureHighlight', {
            "type":"geojson",
            "data": feature.toJSON()
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
      }

      this.setState({
        selectedFeature: feature
      })
    })
  }

  onLayerVisChange(e) {
    console.log(e)
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
            <h4 className="h5">Layers</h4>
            {
              dataLayers.map(layer => {
                return (
                  <div className="form-check" key={'toggleLayer' + layer} >
                    <input className="form-check-input"
                      type="checkbox"
                      data-layer={layer}
                      defaultChecked={true}
                      id={'toggleLayerCheckbox' + layer}
                      onClick={this.onLayerVisChange}/>
                    <label className="form-check-label" htmlFor={'toggleLayerCheckbox' + layer}>
                      {layer}
                    </label>
                  </div>
                )
              })
            }
        </div>

        <FeatureSidebar feature={this.state.selectedFeature} />

        <div className="custom-map-control bottom-right">
          <div>{`Longitude: ${lng.toFixed(2)} Latitude: ${lat.toFixed(2)} Zoom: ${zoom.toFixed(0)}`}</div>
        </div>
        <div ref={el => this.mapContainer = el} className="map" />
      </Fragment>
    );
  }
}

const FeatureSidebar = (props) => {
  if (!props.feature) {
    return null
  }
  return (
    <div className="custom-map-control top-right">
      <h4 className="h5">Selected Feature</h4>
      {JSON.stringify(props.feature.properties)}
    </div>
  )
}

StaticMap.propTypes = {
  map_style: PropTypes.string.isRequired,
  toggleableLayerIds: PropTypes.array,
  clickableLayerAttributes: PropTypes.object
}

export default StaticMap
