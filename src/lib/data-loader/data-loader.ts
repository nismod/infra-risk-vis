export type DataLoaderSubscriber = (loader: DataLoader) => void;

export class DataLoader<T = any> {
  constructor(
    public readonly id: string,
    public readonly layer: string,
    public readonly field: string,
    public readonly fieldParameters: any,
  ) {}

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
    const path = `/api/attributes/${this.field}`;

    const search = new URLSearchParams();
    search.append('layer', this.layer);
    Object.entries(this.fieldParameters).forEach(([k, v]) => search.append(k, v as string));

    const res = await fetch(`${path}?${search}`, {
      method: 'POST',
      body: JSON.stringify(missingIds),
      headers: {
        'Content-Type': 'application/json',
      },
    });

    if (res.status !== 200) {
      return {};
    } else {
      return await res.json();
    }
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
