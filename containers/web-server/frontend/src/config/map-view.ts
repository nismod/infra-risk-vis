interface MapViewConfig {
  initialViewState: {
    latitude: number;
    longitude: number;
    zoom: number;
  };
  viewLimits: {
    minZoom?: number;
    maxZoom?: number;
    minPitch?: number;
    maxPitch?: number;
  };
}

export const mapViewConfig: MapViewConfig = {
  initialViewState: {
    latitude: 20.0,
    longitude: -40.0,
    zoom: 3,
  },
  viewLimits: {
    minZoom: 3,
    maxZoom: 12,
    maxPitch: 0,
  },
};
