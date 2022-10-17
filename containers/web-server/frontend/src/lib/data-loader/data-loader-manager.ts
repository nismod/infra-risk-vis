import stringify from 'json-stable-stringify';

import { FieldSpec } from '@/lib/data-map/view-layers';

import { DataLoader } from './data-loader';

function getLoaderKey(layer: string, fieldSpec: FieldSpec) {
  return `${layer}__${stringify(fieldSpec)}`;
}

export class DataLoaderManager {
  private loaders: { [key: string]: DataLoader } = {};
  private nextLoaderId = 0;

  public getDataLoader(layer: string, fieldSpec: FieldSpec) {
    const loaderKey = getLoaderKey(layer, fieldSpec);
    if (this.loaders[loaderKey] == null) {
      console.log(`Initialising data loader for ${layer} ${stringify(fieldSpec)}`);
      const loader = new DataLoader(this.nextLoaderId.toString(), layer, fieldSpec);
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
