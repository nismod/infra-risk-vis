import _ from 'lodash';
import { selectorFamily } from 'recoil';

import { DataParamGroupConfig, inferDependenciesFromData, inferDomainsFromData } from '@/lib/controls/data-params';

import { HAZARD_DOMAINS_CONFIG } from '@/config/hazards/domains';
import { HazardType } from '@/config/hazards/metadata';

import { rasterSourceDomainsQuery } from './sources';

export const hazardDomainsConfigState = selectorFamily<DataParamGroupConfig, HazardType>({
  key: 'hazardDomainsConfigState',
  get:
    (hazardType: HazardType) =>
    ({ get }) => {
      const { defaults, dependencies, preprocess } = HAZARD_DOMAINS_CONFIG[hazardType];

      const paramNames = Object.keys(defaults);

      // TODO: assumes domains are named the same as hazard types
      const sourceDomains = get(rasterSourceDomainsQuery(hazardType));

      const uniqueDomains = _(sourceDomains)
        // get unique combinations of parameters
        .map((x) => _.pick(x, paramNames))
        .uniqWith(_.isEqual)

        // preprocess all fields according to config
        .map((obj) => _.mapValues(obj, (value, key) => (preprocess[key] ? preprocess[key](value) : value)))
        .value();

      const inferredDomains = inferDomainsFromData(uniqueDomains);
      // sortBy sorts integers correctly, whereas sort stringifies which results in 1, 10, 2, 20 etc
      const sortedDomains = _.mapValues(inferredDomains, (domain) => _.sortBy(domain));

      return {
        paramDomains: sortedDomains,
        paramDefaults: defaults,
        paramDependencies: inferDependenciesFromData(uniqueDomains, dependencies),
      };
    },
});
