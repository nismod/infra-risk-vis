/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { BaseHttpRequest } from './core/BaseHttpRequest';
import type { OpenAPIConfig } from './core/OpenAPI';
import { FetchHttpRequest } from './core/FetchHttpRequest';

import { AttributesService } from './services/AttributesService';
import { ColormapService } from './services/ColormapService';
import { FeaturesService } from './services/FeaturesService';
import { TilesService } from './services/TilesService';

type HttpRequestConstructor = new (config: OpenAPIConfig) => BaseHttpRequest;

export class ApiClient {

    public readonly attributes: AttributesService;
    public readonly colormap: ColormapService;
    public readonly features: FeaturesService;
    public readonly tiles: TilesService;

    public readonly request: BaseHttpRequest;

    constructor(config?: Partial<OpenAPIConfig>, HttpRequest: HttpRequestConstructor = FetchHttpRequest) {
        this.request = new HttpRequest({
            BASE: config?.BASE ?? '',
            VERSION: config?.VERSION ?? '0.1.0',
            WITH_CREDENTIALS: config?.WITH_CREDENTIALS ?? false,
            CREDENTIALS: config?.CREDENTIALS ?? 'include',
            TOKEN: config?.TOKEN,
            USERNAME: config?.USERNAME,
            PASSWORD: config?.PASSWORD,
            HEADERS: config?.HEADERS,
            ENCODE_PATH: config?.ENCODE_PATH,
        });

        this.attributes = new AttributesService(this.request);
        this.colormap = new ColormapService(this.request);
        this.features = new FeaturesService(this.request);
        this.tiles = new TilesService(this.request);
    }
}
