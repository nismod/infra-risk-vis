import { atom } from 'recoil';

import { RegionalExposureVariableType } from '@/config/regional-risk/metadata';

export const regionalExposureVariableState = atom<RegionalExposureVariableType>({
  key: 'regionalExposureVariableState',
  default: 'pop_exposed_seismic_threshold0.1g',
});
