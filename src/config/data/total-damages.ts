import { DataParamConfig } from './types';

interface TotalDamagesParams {
  rcp: string;
  epoch: number;
}

export const totalDamagesConfig: DataParamConfig<TotalDamagesParams> = {
  paramDefaults: {
    rcp: 'baseline',
    epoch: 2010,
  },
  paramDomains: {
    rcp: ['baseline', '4.5', '8.5'],
    epoch: [2010, 2050, 2100],
  },
};
