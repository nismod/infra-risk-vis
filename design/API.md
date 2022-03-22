# API structure

Design document for specifying the structure of data sent from the asset/feature API.

## Asset details endpoint

`GET /assets/{asset_id}`:

```ts
interface AssetDamage {
  hazard: string | 'total';
  rcp: string;
  epoch: number;
  rp: number | 'EA';

  damage_type: 'direct' | 'indirect' | 'combined';
  protection_standard: number | 'undefended';

  min: number;
  mean: number;
  max: number;
}

interface AssetDetailsResponse {
  id: number;
  source_id: string;
  layer: string;

  properties: Record<string, any>;

  damages: AssetDamage[];

  //...
}
```

## Assets list endpoint

`GET /assets?sort_by={sort_field_id}`:

```ts
interface AssetListItem {
  id: number;
  source_id: string;
  layer: string;

  properties: Record<string, any>;
}

type AssetListResponse = AssetListItem[];
```

## Feature extra map data endpoint

`GET /attributes/{attribute_id}/{layer}/{z}/{x}/{y}`:

```ts
type FeatureId = number;

type FeatureAttributesResponse<AttributeType> = {
  z: number;
  x: number;
  y: number;
  values: Record<FeatureId, AttributeType>;
};
```
