/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */

import type { DamageType } from './DamageType';

export type Damage = {
    hazard: string;
    rcp: string;
    epoch: string;
    damage_type: DamageType;
    protection_standard: number;
    min: number;
    mean: number;
    max: number;
};
