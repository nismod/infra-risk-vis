import { selectorFamily } from 'recoil';

import { FeatureOut } from '@/lib/api-client';

import { apiClient } from '@/api-client';

export const apiFeatureQuery = selectorFamily<FeatureOut, number>({
  key: 'apiFeatureQuery',
  get:
    (id: number) =>
    ({ get }) => {
      return apiClient.features.featuresReadFeature({ featureId: id });
    },
});
