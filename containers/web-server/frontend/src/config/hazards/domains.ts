import { toNumber } from 'lodash';

import { HazardType } from './metadata';

export interface HazardParams {
  returnPeriod: number;
  epoch: string;
  rcp: string;
  gcm: string;
}

interface HazardDomainConfig {
  defaults: Record<string, any>;
  dependencies: Record<string, string[]>;
  preprocess: Record<string, (value: any) => any>;
}

const rcpRegex = /\dx\d/;
function preprocessRcp(rcp: string) {
  return rcp.match(rcpRegex) ? rcp.replace('x', '.') : rcp;
}

export const HAZARD_DOMAINS_CONFIG: Record<HazardType, HazardDomainConfig> = {
  fluvial: {
    defaults: {
      rp: 100,
      rcp: 'baseline',
      epoch: 'present',
      gcm: 'MIROC-ESM-CHEM',
    },
    dependencies: {
      rcp: ['epoch'],
      gcm: ['epoch'],
    },
    preprocess: {
      rcp: preprocessRcp,
      rp: toNumber,
    },
  },
  coastal: {
    defaults: {
      rp: 1,
      epoch: '2010',
      rcp: 'baseline',
      gcm: 'None',
    },
    dependencies: {
      rcp: ['epoch'],
    },
    preprocess: {
      rcp: preprocessRcp,
      rp: toNumber,
    },
  },
  cyclone: {
    defaults: {
      rp: 10,
      gcm: 'constant',
      /**
       * epoch and rcp for cyclones are added programmatically upon load
       * adjust custom code for data loading when these fields are added to the backend
       */
      epoch: 'present',
      rcp: 'baseline',
    },
    dependencies: {
      gcm: ['epoch'],
      rcp: ['epoch'],
    },
    preprocess: {
      rp: toNumber,
      rcp: preprocessRcp,
    },
  },
  extreme_heat: {
    defaults: {
      epoch: 'baseline',
      rcp: 'baseline',
      gcm: 'gfdl-esm2m',
    },
    dependencies: {
      rcp: ['epoch'],
    },
    preprocess: {
      rcp: preprocessRcp,
    },
  },
  earthquake: {
    defaults: {
      rp: 475,
      medium: 'soil',
    },
    dependencies: {},
    preprocess: {
      rp: toNumber,
    },
  },
  drought: {
    defaults: {
      epoch: 'baseline',
      rcp: 'baseline',
      gcm: 'gfdl-esm2m',
    },
    dependencies: {
      rcp: ['epoch'],
      gcm: ['epoch'],
    },
    preprocess: {
      rcp: preprocessRcp,
    },
  },
};
