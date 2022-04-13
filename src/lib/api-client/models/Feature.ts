/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */

import type { Damage } from './Damage';

export type Feature = {
    id: number;
    string_id: string;
    layer: string;
    sublayer?: string;
    properties: any;
    damages?: Array<Damage>;
};
