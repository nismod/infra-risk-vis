/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { Feature } from '../models/Feature';

import type { CancelablePromise } from '../core/CancelablePromise';
import type { BaseHttpRequest } from '../core/BaseHttpRequest';

export class FeaturesService {

    constructor(public readonly httpRequest: BaseHttpRequest) {}

    /**
     * Read Feature
     * @param featureId
     * @returns Feature Successful Response
     * @throws ApiError
     */
    public featuresReadFeature(
        featureId: number,
    ): CancelablePromise<Feature> {
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

}