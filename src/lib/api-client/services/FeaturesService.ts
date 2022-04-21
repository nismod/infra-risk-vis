/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { FeatureOut } from '../models/FeatureOut';
import type { Page_FeatureListItemOut_float__ } from '../models/Page_FeatureListItemOut_float__';

import type { CancelablePromise } from '../core/CancelablePromise';
import type { BaseHttpRequest } from '../core/BaseHttpRequest';

export class FeaturesService {

    constructor(public readonly httpRequest: BaseHttpRequest) {}

    /**
     * Read Feature
     * @returns FeatureOut Successful Response
     * @throws ApiError
     */
    public featuresReadFeature({
        featureId,
    }: {
        featureId: number,
    }): CancelablePromise<FeatureOut> {
        return this.httpRequest.request({
            method: 'GET',
            url: '/features/{feature_id}',
            path: {
                'feature_id': featureId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }

    /**
     * Read Sorted Features
     * @returns Page_FeatureListItemOut_float__ Successful Response
     * @throws ApiError
     */
    public featuresReadSortedFeatures({
        field,
        layer,
        fieldParams,
        page = 1,
        size = 50,
    }: {
        field: string,
        layer: string,
        fieldParams: string,
        page?: number,
        size?: number,
    }): CancelablePromise<Page_FeatureListItemOut_float__> {
        return this.httpRequest.request({
            method: 'GET',
            url: '/features/sorted-by/{field}',
            path: {
                'field': field,
            },
            query: {
                'layer': layer,
                'field_params': fieldParams,
                'page': page,
                'size': size,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }

}