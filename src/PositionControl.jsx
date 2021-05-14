import React from 'react';

const PositionControl = (props) => (
  <div className="custom-map-control bottom-right position-control mapboxgl-ctrl-scale">
    Lon: {props.lng.toFixed(2)} Lat: {props.lat.toFixed(2)} Zoom: {props.zoom.toFixed(0)}
  </div>
)

export default PositionControl;
