import { DataParamGroupConfig } from 'lib/controls/data-params';

export interface HazardParams {
  returnPeriod: number;
  epoch: string;
  rcp: string;
  gcm: string;
}

export const HAZARD_DOMAINS: Record<string, DataParamGroupConfig<HazardParams>> = {
  fluvial: {
    paramDomains: {
      returnPeriod: [1, 2, 5, 10, 25, 50, 100, 250, 500, 1000],

      rcp: ['baseline', '2.6', '4.5', '8.5'],
      epoch: ['present', '2030', '2050', '2080'],
      gcm: ["None", "WATCH", "NorESM1-M", "GFDL_ESM2M", "GFDL-ESM2M", "HadGEM2-ES", "IPSL-CM5A-LR", "MIROC-ESM-CHEM"],
    },
    paramDefaults: {
      returnPeriod: 100,
      rcp: 'baseline',
      epoch: 'present',
      gcm: 'MIROC-ESM-CHEM',
    },
    paramDependencies: {
      rcp: ({ epoch }) => {
        if (epoch === 'present') return ['baseline'];
        else if (epoch === '2050' || epoch === '2080') return ['2.6', '4.5', '8.5'];
      },
      gcm: ({ epoch }) => {
        if (epoch === 'present') return ['WATCH'];
        if (epoch !== 'present') return ['GFDL-ESM2M', 'HadGEM2-ES', 'IPSL-CM5A-LR', 'NorESM1-M', 'MIROC-ESM-CHEM'];
      },
    },
  },
  surface: {
    paramDomains: {
      returnPeriod: [20, 50, 100, 200, 500, 1500],

      rcp: ['baseline', '2.6', '4.5', '8.5'],
      epoch: ['2010', '2050', '2080'],
      gcm: ['None'],
    },
    paramDefaults: {
      returnPeriod: 100,

      rcp: 'baseline',
      epoch: '2010',
      gcm: 'None',
    },
    paramDependencies: {
      rcp: ({ epoch }) => {
        if (epoch === '2010') return ['baseline'];
        else if (epoch === '2050' || epoch === '2080') return ['2.6', '4.5', '8.5'];
      },
    },
  },
  coastal: {
    paramDomains: {
      returnPeriod: [1, 2, 5, 10, 50, 100, 250, 500, 1000],
      epoch: ['present', '2030', '2050', '2080'],
      rcp: ['baseline', '2.6', '4.5', '8.5'],

      gcm: ['None'],
    },
    paramDefaults: {
      returnPeriod: 1,
      epoch: '2010',
      rcp: 'baseline',
      gcm: 'None',
    },
    paramDependencies: {
      rcp: ({ epoch }) => {
        if (epoch === '2010') return ['baseline'];
        if (epoch === '2050' || epoch === '2100') return ['2.6', '4.5', '8.5'];
        if (epoch === '2030' || epoch === '2070') return ['4.5', '8.5'];
      },
    },
  },
  cyclone: {
    paramDomains: {
      returnPeriod: [
        10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000,
        6000, 7000, 8000, 9000, 10000,
      ],
      epoch: ['2010', '2050', '2100'],
      rcp: ['baseline', '4.5', '8.5'],
      gcm: ['None'],
    },
    paramDefaults: {
      returnPeriod: 100,
      epoch: '2010',
      rcp: 'baseline',
      gcm: 'None',
    },
    paramDependencies: {
      rcp: ({ epoch }) => {
        if (epoch === '2010') return ['baseline'];
        if (epoch === '2050' || epoch === '2100') return ['4.5', '8.5'];
      },
    },
  },
};
