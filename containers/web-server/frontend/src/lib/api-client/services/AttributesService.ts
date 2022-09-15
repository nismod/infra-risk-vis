/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { CancelablePromise } from '../core/CancelablePromise';
import type { BaseHttpRequest } from '../core/BaseHttpRequest';

export class AttributesService {

    constructor(public readonly httpRequest: BaseHttpRequest) {}

    /**
     * Read Attributes
     * @returns any Successful Response
     * @throws ApiError
     */
    public attributesReadAttributes({
        fieldGroup,
        layer,
        field,
        dimensions,
        parameters,
        requestBody,
    }: {
        fieldGroup: string,
        layer: string,
        field: string,
        dimensions: string,
        parameters: string,
        requestBody: Array<number>,
    }): CancelablePromise<any> {
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
                'parameters': parameters,
            },
            body: requestBody,
            mediaType: 'application/json',
            errors: {
                422: `Validation Error`,
            },
        });
    }

}