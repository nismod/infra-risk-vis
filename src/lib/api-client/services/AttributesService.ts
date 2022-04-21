/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { DamageType } from '../models/DamageType';

import type { CancelablePromise } from '../core/CancelablePromise';
import type { BaseHttpRequest } from '../core/BaseHttpRequest';

export class AttributesService {

    constructor(public readonly httpRequest: BaseHttpRequest) {}

    /**
     * Read Damages
     * @returns number Successful Response
     * @throws ApiError
     */
    public attributesReadDamages({
        layer,
        hazard,
        rcp,
        epoch,
        damageType,
        protectionStandard,
        requestBody,
    }: {
        layer: string,
        hazard: string,
        rcp: string,
        epoch: string,
        damageType: DamageType,
        protectionStandard: number,
        requestBody: Array<number>,
    }): CancelablePromise<Record<string, number>> {
        return this.httpRequest.request({
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