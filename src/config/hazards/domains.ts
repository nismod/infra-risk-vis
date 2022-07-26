import { DataParamGroupConfig } from 'lib/controls/data-params';

export interface HazardParams {
  returnPeriod: number;
  epoch: number;
  rcp: string;
  gcm: string;
}

export const HAZARD_DOMAINS: Record<string, DataParamGroupConfig<HazardParams>> = {
  river: {
    paramDomains: {
      returnPeriod: [2, 5, 10, 25, 50, 100, 250, 500, 1000],
      rcp: ['baseline', '4.5', '8.5'],
      epoch: [2010, 2030, 2050, 2080],
      gcm: ['WATCH', 'GFDL-ESM2M', 'HadGEM2-ES', 'IPSL-CM5A-LR', 'NorESM1-M', 'MIROC-ESM-CHEM'],
    },
    paramDefaults: {
      returnPeriod: 100,
      rcp: 'baseline',
      epoch: 2010,
      gcm: 'WATCH',
    },
    paramDependencies: {
      rcp: ({ epoch }) => {
        if (epoch === 2010) return ['baseline'];
        if (epoch > 2010) return ['4.5', '8.5'];
      },
      gcm: ({ epoch }) => {
        if (epoch === 2010) return ['WATCH'];
        if (epoch > 2010) return ['GFDL-ESM2M', 'HadGEM2-ES', 'IPSL-CM5A-LR', 'NorESM1-M', 'MIROC-ESM-CHEM'];
      },
    },
  },
  coastal: {
    paramDomains: {
      returnPeriod: [2, 5, 10, 25, 50, 100, 250, 500, 1000],
      epoch: [2010, 2030, 2050, 2080],
      rcp: ['baseline', '4.5', '8.5'],
      gcm: ['None'],
    },
    paramDefaults: {
      returnPeriod: 100,
      epoch: 2010,
      rcp: 'baseline',
      gcm: 'None',
    },
    paramDependencies: {
      rcp: ({ epoch }) => {
        if (epoch === 2010) return ['baseline'];
        if (epoch > 2010) return ['4.5', '8.5'];
      },
    },
  },
  cyclone: {
    paramDomains: {
      returnPeriod: [
        10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000
      ],
      epoch: [2010, 2050],
      rcp: ['baseline', '8.5'],
      gcm: ['None', 'CMCC-CM2-VHR4','CNRM-CM6-1-HR', 'EC-Earth3P-HR', 'HADGEM3-GC31-HM'],
    },
    paramDefaults: {
      returnPeriod: 100,
      epoch: 2010,
      rcp: 'baseline',
      gcm: 'None',
    },
    paramDependencies: {
      rcp: ({ epoch }) => {
        if (epoch === 2010) return ['baseline'];
        if (epoch > 2010) return ['8.5'];
      },
      gcm: ({ epoch }) => {
        if (epoch === 2010) return ['None'];
        if (epoch > 2010) return ['CMCC-CM2-VHR4','CNRM-CM6-1-HR', 'EC-Earth3P-HR', 'HADGEM3-GC31-HM'];
      },
    },
  },
};
