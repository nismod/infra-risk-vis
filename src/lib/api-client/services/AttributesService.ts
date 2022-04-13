/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { DamageType } from '../models/DamageType';

import type { CancelablePromise } from '../core/CancelablePromise';
import { OpenAPI } from '../core/OpenAPI';
import { request as __request } from '../core/request';

export class AttributesService {

    /**
     * Read Damages
     * @param layer
     * @param hazard
     * @param rcp
     * @param epoch
     * @param damageType
     * @param protectionStandard
     * @param requestBody
     * @returns number Successful Response
     * @throws ApiError
     */
    public static attributesReadDamages(
        layer: string,
        hazard: string,
        rcp: string,
        epoch: string,
        damageType: DamageType,
        protectionStandard: number,
        requestBody: Array<number>,
    ): CancelablePromise<Record<string, number>> {
        return __request(OpenAPI, {
            method: 'POST',
            url: '/attributes/damages',
            query: {
                'layer': layer,
                'hazard': hazard,
                'rcp': rcp,
                'epoch': epoch,
                'damage_type': damageType,
                'protection_standard': protectionStandard,
            },
            body: requestBody,
            mediaType: 'application/json',
            errors: {
                422: `Validation Error`,
            },
        });
    }

}