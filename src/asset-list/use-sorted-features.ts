import { ApiClient } from 'lib/api-client';
import { FieldSpec } from 'lib/data-map/view-layers';
import _ from 'lodash';
import { useCallback, useEffect, useState } from 'react';

const apiClient = new ApiClient({
  BASE: 'api',
});

interface PageInfo {
  page: number;
  size: number;
  total: number;
}

export const useSortedFeatures = (layer: string, fieldSpec: FieldSpec, page = 1, pageSize = 50) => {
  const [features, setFeatures] = useState([]);
  const [pageInfo, setPageInfo] = useState<PageInfo>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);

  const fetchFeatures = useCallback(async () => {
    setLoading(true);
    setError(null);

    try {
      const { fieldGroup, fieldDimensions, field } = fieldSpec;
      if (fieldGroup !== 'damages') {
        throw new Error('Only damages field is supported');
      }
      const response = await apiClient.features.featuresReadSortedFeatures({
        layer,
        fieldGroup,
        field,
        dimensions: JSON.stringify(fieldDimensions),
        page,
        size: pageSize,
      });
      setFeatures(response.items);

      setPageInfo(_.pick(response, ['page', 'size', 'total']));
    } catch (error) {
      setError(error);
    }

    setLoading(false);
  }, [layer, fieldSpec, page, pageSize]);

  useEffect(() => {
    fetchFeatures();
  }, [fetchFeatures, fieldSpec, page, pageSize]);

  return {
    features,
    pageInfo,
    loading,
    error,
  };
};
