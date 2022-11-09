import _, { toNumber } from 'lodash';
import { selectorFamily } from 'recoil';

import { DataParamGroupConfig, inferDependenciesFromData, inferDomainsFromData } from '@/lib/controls/data-params';

import { HAZARD_DOMAINS_CONFIG } from '@/config/hazards/domains';
import { HazardType } from '@/config/hazards/metadata';

import { rasterSourceDomainsQuery } from './sources';

const rcpRegex = /\dx\d/;
function preprocessRcp(rcp: string) {
  return rcp.match(rcpRegex) ? rcp.replace('x', '.') : rcp;
}

const HAZARD_DOMAIN_PREPROCESSING = {
  rcp: preprocessRcp,
  rp: toNumber,
};

function numbersLast(value) {
  return isNaN(value) ? 0 : 1;
}
/**
 * Last-minute hack for sorting the epoch etc. labels correctly in the UI
 */
const HAZARD_DOMAIN_SORTING_FUNCTIONS = {
  epoch: [numbersLast, _.identity],
  rcp: [numbersLast, _.identity],
};

/**
 * Last-minute hack for adding parameters missing in the backend
 */
const HAZARD_DOMAIN_CUSTOM_MODIFICATIONS = {
  cyclone: (value: any) => ({
    ...value,
    epoch: value.gcm === 'constant' ? 'present' : '2050',

    // match the rcp format for other hazards - it's preprocessed later, anyway
    rcp: value.gcm === 'constant' ? 'baseline' : '8x5',
  }),
};

export const hazardDomainsConfigState = selectorFamily<DataParamGroupConfig, HazardType>({
  key: 'hazardDomainsConfigState',
  get:
    (hazardType: HazardType) =>
    ({ get }) => {
      const { defaults, dependencies } = HAZARD_DOMAINS_CONFIG[hazardType];
      const preprocess = HAZARD_DOMAIN_PREPROCESSING;

      const paramNames = Object.keys(defaults);

      const customTransform = HAZARD_DOMAIN_CUSTOM_MODIFICATIONS[hazardType] ?? _.identity;

      // TODO: assumes domains are named the same as hazard types
      const sourceDomains = get(rasterSourceDomainsQuery(hazardType)).map(customTransform);

      const uniqueDomains = _(sourceDomains)
        // get unique combinations of parameters
        .map((x) => _.pick(x, paramNames))
        .uniqWith(_.isEqual)

        // preprocess all fields according to config
        .map((obj) => _.mapValues(obj, (value, key) => (preprocess[key] ? preprocess[key](value) : value)))
        .value();

      const inferredDomains = inferDomainsFromData(uniqueDomains);
      // sortBy sorts integers correctly, whereas sort stringifies which results in 1, 10, 2, 20 etc
      const sortedDomains = _.mapValues(inferredDomains, (domain, domainName) =>
        _.sortBy(domain, HAZARD_DOMAIN_SORTING_FUNCTIONS[domainName]),
      );
      return {
        paramDomains: sortedDomains,
        paramDefaults: defaults,
        paramDependencies: inferDependenciesFromData(uniqueDomains, dependencies, HAZARD_DOMAIN_SORTING_FUNCTIONS),
      };
    },
});
