/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
export { ApiClient } from './ApiClient';

export { ApiError } from './core/ApiError';
export { BaseHttpRequest } from './core/BaseHttpRequest';
export { CancelablePromise, CancelError } from './core/CancelablePromise';
export { OpenAPI } from './core/OpenAPI';
export type { OpenAPIConfig } from './core/OpenAPI';

export type { Damage } from './models/Damage';
export { DamageType } from './models/DamageType';
export type { Feature } from './models/Feature';
export type { HTTPValidationError } from './models/HTTPValidationError';
export type { ValidationError } from './models/ValidationError';

export { AttributesService } from './services/AttributesService';
export { FeaturesService } from './services/FeaturesService';
