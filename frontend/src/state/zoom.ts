import { mapViewConfig } from 'config/map-view';
import { atom } from 'recoil';

export const zoomState = atom({
  key: 'zoom',
  default: mapViewConfig.initialViewState.zoom,
});
