import { DataParamGroupConfig } from 'lib/controls/data-params';

export interface HazardParams {
  returnPeriod: number;
  epoch: number;
  rcp: string;
}

export const HAZARD_DOMAINS: Record<string, DataParamGroupConfig<HazardParams>> = {
  river: {
    paramDomains: {
      returnPeriod: [2, 5, 10, 25, 50, 100, 250, 500, 1000],
      rcp: ['baseline', '4.5', '8.5'],
      epoch: [2010, 2030, 2050, 2080],
    },
    paramDefaults: {
      returnPeriod: 100,

      rcp: 'baseline',
      epoch: 2010,
    },
    paramDependencies: {
      rcp: ({ epoch }) => {
        if (epoch === 2010) return ['baseline'];
        else return ['4.5', '8.5'];
      },
    },
  },
  coastal: {
    paramDomains: {
      returnPeriod: [2, 5, 10, 25, 50, 100, 250, 500, 1000],
      rcp: ['baseline', '4.5', '8.5'],
      epoch: [2010, 2030, 2050, 2080],
    },
    paramDefaults: {
      returnPeriod: 100,
      epoch: 2010,
      rcp: 'baseline',
    },
    paramDependencies: {
      rcp: ({ epoch }) => {
        if (epoch === 2010) return ['baseline'];
        else return ['4.5', '8.5'];
      },
    },
  },
};
