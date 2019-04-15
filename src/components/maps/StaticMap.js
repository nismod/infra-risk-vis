import React from 'react'
import PropTypes from 'prop-types'

import mapboxgl from 'mapbox-gl'

class StaticMap extends React.Component {
  map

  constructor(props) {
    super(props);
    this.state = {
      lng: -56,
      lat: -30,
      zoom: 3
    };
    this.map = undefined

    this.onLayerVisClick = this.onLayerVisClick.bind(this);
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

  render() {
    const { lng, lat, zoom } = this.state
    const { toggleableLayerIds } = this.props

    return (
      <div>
        <div className="absolute top left mt12 ml12 bg-darken75 color-white z1 py6 px12 txt-s txt-bold">
          <div className="mt6 mb12">
            <div>Layers</div>
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
  toggleableLayerIds: PropTypes.array
}

export default StaticMap