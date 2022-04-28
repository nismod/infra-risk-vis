import { DataParamGroupConfig, inferDependenciesFromData, inferDomainsFromData } from 'lib/controls/data-params';

import adaptationOptions from './adaptation-options.json';

export interface AdaptationOptionParams {
  sector: string;
  subsector: string;
  hazard: string;
  rcp: string;
  adaptation_name: string;
  adaptation_protection_level: number;
}

export const adaptationDomainsConfig: DataParamGroupConfig<AdaptationOptionParams> = {
  paramDefaults: {
    sector: 'power',
    subsector: 'transmission',
    hazard: 'flooding',
    rcp: '2.6',
    adaptation_name: 'Building protective wall',
    adaptation_protection_level: 1,
  },
  paramDomains: inferDomainsFromData<AdaptationOptionParams>(adaptationOptions),
  paramDependencies: inferDependenciesFromData(adaptationOptions, {
    sector: [],
    subsector: ['sector'],
    hazard: ['sector', 'subsector'],
    rcp: ['sector', 'subsector'],
    adaptation_name: ['sector', 'subsector', 'hazard'],
    adaptation_protection_level: ['sector', 'subsector', 'hazard', 'adaptation_name'],
  }),
};
