import { selector, selectorFamily } from 'recoil';

import { apiClient } from '@/api-client';

export const rasterAllSourcesQuery = selector({
  key: 'rasterAllSourcesQuery',
  get: async () => {
    return await apiClient.tiles.tilesGetAllTileSourceMeta();
  },
});

export const rasterSourceByDomainQuery = selectorFamily({
  key: 'rasterSourceByDomainQuery',
  get:
    (domain: string) =>
    ({ get }) => {
      const sources = get(rasterAllSourcesQuery);

      const sourcesWithDomain = sources.filter((x) => x.domain === domain);

      if (sourcesWithDomain.length > 1) {
        throw new Error(`More than one raster source with domain: ${domain}`);
      }
      if (sourcesWithDomain.length === 0) {
        throw new Error(`No raster sources found with domain: ${domain}`);
      }

      return sourcesWithDomain[0];
    },
});

export const rasterSourceDomainsQuery = selectorFamily({
  key: 'rasterSourceDomainsQuery',
  get:
    (domain: string) =>
    async ({ get }) => {
      const { id } = get(rasterSourceByDomainQuery(domain));

      const { domains } = await apiClient.tiles.tilesGetTileSourceDomains({ sourceId: id });

      return domains;
    },
});
