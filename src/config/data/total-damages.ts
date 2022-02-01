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
  paramDependencies: {
    rcp: ({ epoch }) => {
      if (epoch === 2010) return ['baseline'];
      if (epoch === 2050 || epoch === 2100) return ['4.5', '8.5'];
    },
  },
};
