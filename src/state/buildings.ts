import { atom } from 'recoil';

export const buildingsStyleState = atom({
  key: 'buildingsStyleState',
  default: 'type',
});

export const buildingSelectionState = atom({
  key: 'buildingSelectionState',
  default: {
    buildings_commercial: true,
    buildings_residential: true,
    buildings_institutional: true,
    buildings_mixed: true,
    buildings_industrial: true,
    buildings_recreation: true,
    buildings_other: true,
    buildings_resort: true,
  },
});
