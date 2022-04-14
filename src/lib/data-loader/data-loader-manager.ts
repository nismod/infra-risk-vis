import { FieldSpec } from 'lib/data-map/view-layers';
import { DataLoader } from './data-loader';

function getLoaderKey(layer: string, fieldSpec: FieldSpec) {
  return `${layer}__${fieldSpec.field}__${Object.entries(fieldSpec.fieldParams)
    .map(([k, v]) => `${k}_${v}`)
    .join('__')}`;
}

export class DataLoaderManager {
  private loaders: { [key: string]: DataLoader } = {};
  private nextLoaderId = 0;

  public getDataLoader(layer: string, fieldSpec: FieldSpec) {
    const loaderKey = getLoaderKey(layer, fieldSpec);
    if (this.loaders[loaderKey] == null) {
      const loader = new DataLoader(this.nextLoaderId.toString(), layer, fieldSpec.field, fieldSpec.fieldParams);
      this.nextLoaderId += 1;
      this.loaders[loaderKey] = loader;
    }
    return this.loaders[loaderKey];
  }

  public clearDataLoader(layer: string, fieldSpec: FieldSpec) {
    const loaderKey = getLoaderKey(layer, fieldSpec);

    if (this.loaders[loaderKey] != null) {
      this.loaders[loaderKey].destroy();
      delete this.loaders[loaderKey];
    }
  }
}

export const dataLoaderManager = new DataLoaderManager();
