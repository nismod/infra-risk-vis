/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { CancelablePromise } from '../core/CancelablePromise';
import type { BaseHttpRequest } from '../core/BaseHttpRequest';

export class AttributesService {

    constructor(public readonly httpRequest: BaseHttpRequest) {}

    /**
     * Read Attributes
     * @returns number Successful Response
     * @throws ApiError
     */
    public attributesReadAttributes({
        fieldGroup,
        layer,
        field,
        dimensions,
        requestBody,
    }: {
        fieldGroup: string,
        layer: string,
        field: string,
        dimensions: string,
        requestBody: Array<number>,
    }): CancelablePromise<Record<string, number>> {
        return this.httpRequest.request({
            method: 'POST',
            url: '/attributes/{field_group}',
            path: {
                'field_group': fieldGroup,
            },
            query: {
                'layer': layer,
                'field': field,
                'dimensions': dimensions,
            },
            body: requestBody,
            mediaType: 'application/json',
            errors: {
                422: `Validation Error`,
            },
        });
    }

}