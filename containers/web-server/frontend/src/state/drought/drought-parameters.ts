import { atom } from 'recoil';

import { DroughtOptionsVariableType, DroughtRiskVariableType } from '@/config/drought/metadata';

export const droughtRcpParamState = atom<string>({
  key: 'droughtRcpParamState',
  default: '2.6',
});

export const droughtShowRiskState = atom<boolean>({
  key: 'droughtShowRiskState',
  default: true,
});

export const droughtRiskVariableState = atom<DroughtRiskVariableType>({
  key: 'droughtRiskVariableState',
  default: 'mean_monthly_water_stress_',
});

export const droughtShowOptionsState = atom<boolean>({
  key: 'droughtShowOptionsState',
  default: true,
});

export const droughtOptionsVariableState = atom<DroughtOptionsVariableType>({
  key: 'droughtOptionsVariableState',
  default: 'cost_jmd',
});
