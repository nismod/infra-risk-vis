import { DataParamGroupConfig, inferDependenciesFromData, inferDomainsFromData } from '@/lib/controls/data-params';

import adaptationOptions from './adaptation-options.json';

export interface AdaptationOptionParams {
  sector: string;
  subsector: string;
  asset_type: string;
  hazard: string;
  rcp: string;
  adaptation_name: string;
  adaptation_protection_level: number;
}

export const adaptationDomainsConfig: DataParamGroupConfig<AdaptationOptionParams> = {
  paramDefaults: {
    sector: 'power',
    subsector: 'transmission',
    asset_type: 'pole',
    hazard: 'flooding',
    rcp: '2.6',
    adaptation_name: 'Building protective wall',
    adaptation_protection_level: 1,
  },
  paramDomains: inferDomainsFromData<AdaptationOptionParams>(adaptationOptions),
  paramDependencies: inferDependenciesFromData(adaptationOptions, {
    sector: [],
    subsector: ['sector'],
    asset_type: ['sector', 'subsector'],
    hazard: ['sector', 'subsector', 'asset_type'],
    rcp: ['sector', 'subsector', 'asset_type'],
    adaptation_name: ['sector', 'subsector', 'asset_type', 'hazard'],
    adaptation_protection_level: ['sector', 'subsector', 'asset_type', 'hazard', 'adaptation_name'],
  }),
};
