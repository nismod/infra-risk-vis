/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { Feature } from '../models/Feature';

import type { CancelablePromise } from '../core/CancelablePromise';
import { OpenAPI } from '../core/OpenAPI';
import { request as __request } from '../core/request';

export class FeaturesService {

    /**
     * Read Feature
     * @param featureId
     * @returns Feature Successful Response
     * @throws ApiError
     */
    public static featuresReadFeature(
        featureId: number,
    ): CancelablePromise<Feature> {
        return __request(OpenAPI, {
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

}