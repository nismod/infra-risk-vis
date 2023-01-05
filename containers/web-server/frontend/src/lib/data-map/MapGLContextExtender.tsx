import _ from 'lodash';
import { FC, useContext, useMemo } from 'react';
import { MapContext } from 'react-map-gl';

export interface MapGLContextExtenderProps {
  viewLimits: {
    minZoom?: number;
    maxZoom?: number;
    minPitch?: number;
    maxPitch?: number;
    minBearing?: number;
    maxBearing?: number;
  };
}

/**
 * Needed to add the missing view state limits (maxZoom etc) to the react-map-gl MapContext.
 *
 * Deck.gl doesn't forward these properties to the context it creates,
 * so react-map-gl components like NavigationControl end up resetting the viewLimits to defaults.
 */
export const MapGLContextExtender: FC<MapGLContextExtenderProps> = ({ viewLimits, children }) => {
  const baseContext = useContext(MapContext);

  // need to make a plain object out of viewport, and then assign the missing fields
  const plainViewport = Object.fromEntries(Object.entries(baseContext.viewport ?? {}));
  const viewport: any = Object.assign(plainViewport, viewLimits);

  const extendedContext = useMemo(
    () => ({ ..._.omit(baseContext, 'viewport'), viewport }),
    [baseContext, viewport],
  );
  return <MapContext.Provider value={extendedContext}>{children}</MapContext.Provider>;
};
