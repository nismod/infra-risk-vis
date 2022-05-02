/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */

import type { Adaptation } from './Adaptation';
import type { ExpectedDamage } from './ExpectedDamage';
import type { NPVDamage } from './NPVDamage';
import type { ReturnPeriodDamage } from './ReturnPeriodDamage';

export type FeatureOut = {
    id: number;
    string_id: string;
    layer: string;
    sublayer?: string;
    properties: any;
    damages_expected?: Array<ExpectedDamage>;
    damages_return_period?: Array<ReturnPeriodDamage>;
    damages_npv?: Array<NPVDamage>;
    adaptation?: Array<Adaptation>;
};
