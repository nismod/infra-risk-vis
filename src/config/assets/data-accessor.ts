import { DataLoader } from 'lib/data-loader/data-loader';
import { DataManager } from 'lib/data-map/view-layers';
import { extraProperty, featureProperty } from 'lib/deck/props/data-source';

export interface FieldSpec {
  field: string;
  fieldParams?: any;
}

function getLoaderKey(layer: string, fieldSpec: FieldSpec) {
  return `${layer}__${fieldSpec.field}__${Object.entries(fieldSpec.fieldParams)
    .map(([k, v]) => `${k}_${v}`)
    .join('__')}`;
}

export class AssetDataManager implements DataManager {
  private loaders: { [key: string]: DataLoader } = {};

  public getDataAccessor(layer: string, fieldSpec: FieldSpec) {
    if (typeof fieldSpec === 'string' || typeof fieldSpec === 'function') {
      return featureProperty(fieldSpec);
    } else {
      const dataLoader = this.getDataLoader(layer, fieldSpec);
      return extraProperty(dataLoader);
    }
  }

  public getDataLoader(layer: string, fieldSpec: FieldSpec) {
    const loaderKey = getLoaderKey(layer, fieldSpec);
    if (this.loaders[loaderKey] == null) {
      const loader = new DataLoader(layer, fieldSpec.field, fieldSpec.fieldParams);
      this.loaders[loaderKey] = loader;
      // loader.subscribe();
    }
    return this.loaders[loaderKey];
  }
}

export const assetDataManager = new AssetDataManager();
