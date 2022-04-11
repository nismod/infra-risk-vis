// function setAtPath(obj: any, path: any[], value: any) {
//   for (const prop of path.slice(0, path.length - 1)) {
//     obj[prop] ??= {};
//     obj = obj[prop];
//   }
//   obj[path[path.length - 1]] = value;
// }

export type DataLoaderSubscriber = (loader: DataLoader) => void;

export class DataLoader<T = any> {
  constructor(layer: string, field: string, fieldParameters: any) {
    this.layer = layer;
    this.field = field;
    this.fieldParameters = fieldParameters;
  }

  public readonly layer: string;
  public readonly field: string;
  public readonly fieldParameters: any;

  private _updateTrigger: number = 1;

  public get updateTrigger() {
    return this._updateTrigger;
  }

  private data: Map<number, T> = new Map();
  // private loadedTilesZXY: Record<number, Record<number, Record<number, boolean>>> = {};
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

  // private async requestTileData(x, y, z): Promise<Record<number, T>> {
  //   // MOCK DATA
  //   if (
  //     this.layer === 'elec_edges_high' &&
  //     this.field === 'damages' &&
  //     this.fieldParameters?.damageType === 'direct' &&
  //     this.fieldParameters?.hazard === 'cyclone' &&
  //     this.fieldParameters?.rcp === '8.5' &&
  //     this.fieldParameters?.epoch === '2100'
  //   ) {
  //     return new Promise((resolve, reject) => {
  //       setTimeout(() => {
  //         resolve(mockElecHighDamages[z]?.[x]?.[y] ?? {});
  //       }, 1000 + Math.random() * 3000);
  //     });
  //   } else {
  //     return Promise.resolve({});
  //   }
  //   // const res = await fetch(`/api/${z}/${x}/${y}`); // TODO
  //   // return res.json();
  // }

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

  // async loadTileData(x: number, y: number, z: number) {
  //   const loaded = this.loadedTilesZXY[z]?.[x]?.[y];
  //   if (loaded) return;

  //   const loadedData = await this.requestTileData(x, y, z);
  //   setAtPath(this.loadedTilesZXY, [z, x, y], true);

  //   let newData = false;
  //   for (const [key, value] of Object.entries(loadedData)) {
  //     if (this.data[key] == null) {
  //       newData = true;
  //       this.data[key] = value;
  //     }
  //   }

  //   if (newData) {
  //     this._updateTrigger += 1;
  //     this.subscribers?.forEach((subscriber) => subscriber(this));
  //     console.log('Updated extra data', this._updateTrigger);
  //   }
  // }
}
