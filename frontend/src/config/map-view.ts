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
    latitude: 18.14,
    longitude: -77.28,
    zoom: 8,
  },
  viewLimits: {
    minZoom: 3,
    maxZoom: 16,
    maxPitch: 0,
  },
};
