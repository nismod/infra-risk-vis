import React from 'react';
import PropTypes from 'prop-types';

const PositionControl = ({ lat, lng, zoom }) => (
  <div className="custom-map-control bottom-right position-control mapboxgl-ctrl-scale">
    Lon: {lng.toFixed(2)} Lat: {lat.toFixed(2)} Zoom: {zoom.toFixed(0)}
  </div>
);

PositionControl.propTypes = {
  lat: PropTypes.number.isRequired,
  lng: PropTypes.number.isRequired,
  zoom: PropTypes.number.isRequired,
};

export default PositionControl;
