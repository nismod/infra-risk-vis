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
  private missingIds: Set<number> = new Set();

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
  }

  async loadMissingData() {
    if (this.missingIds.size === 0) return;

    const tempMissingIds = Array.from(this.missingIds);
    const loadedData = await this.requestMissingData(tempMissingIds);

    this.updateData(loadedData);
  }

  async loadDataForIds(ids: number[]) {
    const tempMissingIds = ids.filter((id) => this.data.get(id) === undefined);
    if (tempMissingIds.length === 0) return;

    const loadedData = await this.requestMissingData(tempMissingIds);
    this.updateData(loadedData);
  }

  private async requestMissingData(missingIds: number[]): Promise<Record<string, T>> {
    const { fieldGroup, field, fieldDimensions } = this.fieldSpec;
    console.log(`Requesting missing data`, JSON.stringify(fieldDimensions), missingIds);
    return await apiClient.attributes.attributesReadAttributes({
      layer: this.layer,
      fieldGroup,
      field,
      dimensions: JSON.stringify(fieldDimensions),
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
    }

    if (newData) {
      this._updateTrigger += 1;
      this.subscribers?.forEach((subscriber) => subscriber(this));
    }
  }
}
