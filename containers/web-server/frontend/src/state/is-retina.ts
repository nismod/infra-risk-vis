import { atom } from 'recoil';

function checkIsRetina() {
  // taken from Leaflet source: https://github.com/Leaflet/Leaflet/blob/ee71642691c2c71605bacff69456760cfbc80a2a/src/core/Browser.js#L119
  return (
    (window.devicePixelRatio ||
      (window.screen as any).deviceXDPI / (window.screen as any).logicalXDPI) > 1
  );
}

export const isRetinaState = atom({
  key: 'isRetinaState',
  default: checkIsRetina(),
});
