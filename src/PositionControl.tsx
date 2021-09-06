import React, { FC } from 'react';

interface PositionControlProps {
  lat: number;
  lng: number;
  zoom: number;
}

const PositionControl: FC<PositionControlProps> = ({ lat, lng, zoom }) => (
  <div className="custom-map-control bottom-right position-control mapboxgl-ctrl-scale">
    Lon: {lng.toFixed(2)} Lat: {lat.toFixed(2)} Zoom: {zoom.toFixed(0)}
  </div>
);

export default PositionControl;
