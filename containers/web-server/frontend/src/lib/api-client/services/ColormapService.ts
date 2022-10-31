/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { ColorMap } from '../models/ColorMap';

import type { CancelablePromise } from '../core/CancelablePromise';
import type { BaseHttpRequest } from '../core/BaseHttpRequest';

export class ColormapService {

    constructor(public readonly httpRequest: BaseHttpRequest) {}

    /**
     * Get Colormap
     * Retrieve colormap
     * @returns ColorMap Successful Response
     * @throws ApiError
     */
    public colormapGetColormap({
        colormap,
        stretchRange,
        numValues = 255,
    }: {
        colormap: string,
        stretchRange: string,
        numValues?: number,
    }): CancelablePromise<ColorMap> {
        return this.httpRequest.request({
            method: 'GET',
            url: '/colormap',
            query: {
                'colormap': colormap,
                'stretch_range': stretchRange,
                'num_values': numValues,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }

}