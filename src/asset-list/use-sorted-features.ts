import { ApiClient } from 'lib/api-client';
import { FieldSpec } from 'lib/data-map/view-layers';
import _ from 'lodash';
import { useEffect, useState } from 'react';

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

  const fetchFeatures = async () => {
    setLoading(true);
    setError(null);

    try {
      const { field, fieldParams } = fieldSpec;
      if (field !== 'damages') {
        throw new Error('Only damages field is supported');
      }
      const response = await apiClient.features.featuresReadSortedFeatures({
        layer,
        field,
        fieldParams: JSON.stringify(fieldParams),
        page,
        size: pageSize,
      });
      setFeatures(response.items);

      setPageInfo(_.pick(response, ['page', 'size', 'total']));
    } catch (error) {
      setError(error);
    }

    setLoading(false);
  };

  useEffect(() => {
    fetchFeatures();
  }, [fieldSpec, page, pageSize]);

  return {
    features,
    pageInfo,
    loading,
    error,
  };
};
