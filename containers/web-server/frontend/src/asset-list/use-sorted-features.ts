import { parseSync } from '@loaders.gl/core';
import { WKTLoader } from '@loaders.gl/wkt';
import bbox from '@turf/bbox';
import _ from 'lodash';
import { useCallback, useEffect, useState } from 'react';

import { FeatureListItemOut_float_ } from '@/lib/api-client';
import { BoundingBox } from '@/lib/bounding-box';
import { FieldSpec } from '@/lib/data-map/view-layers';

import { apiClient } from '@/api-client';

export interface PageInfo {
  page: number;
  size: number;
  total: number;
}

export interface LayerSpec {
  layer?: string;
  sector?: string;
  subsector?: string;
  assetType?: string;
}
export type ListFeature = Omit<FeatureListItemOut_float_, 'bbox_wkt'> & {
  bbox: BoundingBox;
};

function processFeature(f: FeatureListItemOut_float_): ListFeature {
  const originalBboxGeom = parseSync(f.bbox_wkt, WKTLoader);
  const processedBbox: BoundingBox = bbox(originalBboxGeom) as BoundingBox;

  return {
    ...f,
    bbox: processedBbox,
  };
}

export const useSortedFeatures = (
  layerSpec: LayerSpec,
  fieldSpec: FieldSpec,
  page = 1,
  pageSize = 50,
) => {
  const [features, setFeatures] = useState([]);
  const [pageInfo, setPageInfo] = useState<PageInfo>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);

  const fetchFeatures = useCallback(async () => {
    setLoading(true);
    setError(null);

    try {
      const { fieldGroup, fieldDimensions, field, fieldParams } = fieldSpec;
      // if (fieldGroup !== 'damages') {
      //   throw new Error('Only damages field is supported');
      // }
      const response = await apiClient.features.featuresReadSortedFeatures({
        ...layerSpec,
        fieldGroup,
        field,
        dimensions: JSON.stringify(fieldDimensions),
        parameters: JSON.stringify(fieldParams),
        page,
        size: pageSize,
      });
      const features = (response.items as FeatureListItemOut_float_[]).map(processFeature);
      setFeatures(features);

      setPageInfo(_.pick(response, ['page', 'size', 'total']));
    } catch (error) {
      setError(error);
    }

    setLoading(false);
  }, [layerSpec, fieldSpec, page, pageSize]);

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
