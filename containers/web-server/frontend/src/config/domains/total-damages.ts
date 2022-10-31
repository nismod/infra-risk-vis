import { DataParamGroupConfig } from '@/lib/controls/data-params';

interface TotalDamagesParams {
  rcp: string;
  epoch: number;
}

export const totalDamagesConfig: DataParamGroupConfig<TotalDamagesParams> = {
  paramDefaults: {
    rcp: 'baseline',
    epoch: 2010,
  },
  paramDomains: {
    rcp: ['baseline', '4.5', '8.5'],
    epoch: [2010, 2050],
  },
  paramDependencies: {
    rcp: ({ epoch }) => {
      if (epoch === 2010) return ['baseline'];
      if (epoch === 2050) return ['4.5', '8.5'];
    },
  },
};
