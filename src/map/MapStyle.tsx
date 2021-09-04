import { useContext, useEffect } from 'react';
import { MapContext } from 'react-map-gl';

export const MapStyle = ({ mapStyle }) => {
  const { map } = useContext(MapContext);

  useEffect(() => {
    map.setStyle(mapStyle, { diff: true });
  }, [map, mapStyle]);

  return <></>;
};
