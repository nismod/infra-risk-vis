import { ApiClient } from 'lib/api-client';
import { FieldSpec } from 'lib/data-map/view-layers';

export type DataLoaderSubscriber = (loader: DataLoader) => void;

const apiClient = new ApiClient({
  BASE: '/api',
});
export class DataLoader<T = any> {
  constructor(public readonly id: string, public readonly layer: string, public readonly fieldSpec: FieldSpec) {}

  private _updateTrigger: number = 1;

  public get updateTrigger() {
    return this._updateTrigger;
  }

  private data: Map<number, T> = new Map();

  // feature IDs that have not been loaded yet
  private missingIds: Set<number> = new Set();

  // feature IDs that are currently being loaded
  private loadingIds: Set<number> = new Set();

  private subscribers: DataLoaderSubscriber[];

  getData(id: number) {
    const data = this.data.get(id);

    if (data === undefined) {
      this.missingIds.add(id);
    }

    return data;
  }

  subscribe(callback: DataLoaderSubscriber) {
    this.subscribers ??= [];
    this.subscribers.push(callback);
  }

  unsubscribe(callback: DataLoaderSubscriber) {
    this.subscribers = this.subscribers?.filter((subscriber) => subscriber !== callback);
  }

  destroy() {
    this.subscribers = [];
    this.data.clear();
    this.missingIds.clear();
    this.loadingIds.clear();
  }

  async loadMissingData() {
    if (this.missingIds.size === 0) return;

    const tempMissingIds = Array.from(this.missingIds).filter((id) => !this.loadingIds.has(id));

    if (tempMissingIds.length === 0) return;

    const loadedData = await this.requestMissingData(tempMissingIds);

    this.updateData(loadedData);
  }

  async loadDataForIds(ids: number[]) {
    const tempMissingIds = ids.filter((id) => this.data.get(id) === undefined && !this.loadingIds.has(id));
    if (tempMissingIds.length === 0) return;

    const loadedData = await this.requestMissingData(tempMissingIds);
    this.updateData(loadedData);
  }

  private async requestMissingData(requestedIds: number[]): Promise<Record<string, T>> {
    const { fieldGroup, field, fieldDimensions, fieldParams } = this.fieldSpec;
    const missingIds = requestedIds.filter((id) => !this.loadingIds.has(id));

    if (missingIds.length === 0) return {};

    missingIds.forEach((id) => this.loadingIds.add(id));

    return await apiClient.attributes.attributesReadAttributes({
      layer: this.layer,
      fieldGroup,
      field,
      dimensions: JSON.stringify(fieldDimensions),
      parameters: JSON.stringify(fieldParams),
      requestBody: missingIds,
    });
  }

  private updateData(loadedData: Record<string, T>) {
    let newData = false;
    for (const [key, value] of Object.entries(loadedData)) {
      const numKey = parseInt(key, 10);
      if (this.data.get(numKey) === undefined) {
        newData = true;
        this.data.set(numKey, value);
      }

      this.missingIds.delete(numKey);
      this.loadingIds.delete(numKey);
    }

    if (newData) {
      this._updateTrigger += 1;
      this.subscribers?.forEach((subscriber) => subscriber(this));
    }
  }
}
