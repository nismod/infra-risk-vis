/* istanbul ignore file */
/* tslint:disable */
/* eslint-disable */
import type { TileSourceDomains } from '../models/TileSourceDomains';
import type { TileSourceMeta } from '../models/TileSourceMeta';

import type { CancelablePromise } from '../core/CancelablePromise';
import type { BaseHttpRequest } from '../core/BaseHttpRequest';

export class TilesService {

    constructor(public readonly httpRequest: BaseHttpRequest) {}

    /**
     * Get All Tile Source Meta
     * Retrieve metadata about all the tile sources available
     * @returns TileSourceMeta Successful Response
     * @throws ApiError
     */
    public tilesGetAllTileSourceMeta(): CancelablePromise<Array<TileSourceMeta>> {
        return this.httpRequest.request({
            method: 'GET',
            url: '/tiles/sources',
        });
    }

    /**
     * Insert Source Meta
     * Ingest Tile Source Meta
     * @returns any Successful Response
     * @throws ApiError
     */
    public tilesInsertSourceMeta({
        xToken,
        requestBody,
    }: {
        xToken: string,
        requestBody: TileSourceMeta,
    }): CancelablePromise<any> {
        return this.httpRequest.request({
            method: 'POST',
            url: '/tiles/sources',
            headers: {
                'x-token': xToken,
            },
            body: requestBody,
            mediaType: 'application/json',
            errors: {
                422: `Validation Error`,
            },
        });
    }

    /**
     * Get Tile Source Meta
     * Retrieve metadata about a single tile source
     * @returns TileSourceMeta Successful Response
     * @throws ApiError
     */
    public tilesGetTileSourceMeta({
        sourceId,
    }: {
        sourceId: number,
    }): CancelablePromise<TileSourceMeta> {
        return this.httpRequest.request({
            method: 'GET',
            url: '/tiles/sources/{source_id}',
            path: {
                'source_id': sourceId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }

    /**
     * Delete Source Meta
     * Delete Tile Source Meta
     * @returns any Successful Response
     * @throws ApiError
     */
    public tilesDeleteSourceMeta({
        sourceId,
        xToken,
    }: {
        sourceId: number,
        xToken: string,
    }): CancelablePromise<any> {
        return this.httpRequest.request({
            method: 'DELETE',
            url: '/tiles/sources/{source_id}',
            path: {
                'source_id': sourceId,
            },
            headers: {
                'x-token': xToken,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }

    /**
     * Get Tile Source Domains
     * Retrieve all combinations available for the source domain
     * @returns TileSourceDomains Successful Response
     * @throws ApiError
     */
    public tilesGetTileSourceDomains({
        sourceId,
    }: {
        sourceId: number,
    }): CancelablePromise<TileSourceDomains> {
        return this.httpRequest.request({
            method: 'GET',
            url: '/tiles/sources/{source_id}/domains',
            path: {
                'source_id': sourceId,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }

    /**
     * Get Tile
     * This route does not work in a FastAPI thread pool environment (i.e. when not async)
     * @returns any Successful Response
     * @throws ApiError
     */
    public tilesGetTile({
        keys,
        tileZ,
        tileX,
        tileY,
        colormap,
        stretchRange,
    }: {
        keys: string,
        tileZ: number,
        tileX: number,
        tileY: number,
        colormap?: string,
        stretchRange?: string,
    }): CancelablePromise<any> {
        return this.httpRequest.request({
            method: 'GET',
            url: '/tiles/{keys}/{tile_z}/{tile_x}/{tile_y}.png',
            path: {
                'keys': keys,
                'tile_z': tileZ,
                'tile_x': tileX,
                'tile_y': tileY,
            },
            query: {
                'colormap': colormap,
                'stretch_range': stretchRange,
            },
            errors: {
                422: `Validation Error`,
            },
        });
    }

}