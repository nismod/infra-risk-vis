import React from 'react';

const PositionControl = (props) => (
  <div className="custom-map-control bottom-right">
    Longitude: {props.lng.toFixed(2)}
    Latitude: {props.lat.toFixed(2)}
    Zoom: {props.zoom.toFixed(0)}
  </div>
)

export default PositionControl;
