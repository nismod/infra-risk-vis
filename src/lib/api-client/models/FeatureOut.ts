/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */

import type { DamageOut } from './DamageOut';

export type FeatureOut = {
    id: number;
    string_id: string;
    layer: string;
    sublayer?: string;
    properties: any;
    damages?: Array<DamageOut>;
};
