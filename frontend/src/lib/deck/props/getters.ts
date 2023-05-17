import { DataLoader } from 'lib/data-loader/data-loader';

export type Accessor<To, From = any> = ((x: From) => To) & { updateTriggers?: any[]; dataLoader?: DataLoader };
export type Getter<T> = T | Accessor<T>;

export function mergeTriggers(...accessors: Accessor<any>[]) {
  const res = [];
  for (const acc of accessors) {
    for (const elem of acc.updateTriggers ?? []) {
      res.push(elem);
    }
  }
  return res;
}

export function withTriggers(fn: Accessor<any>, triggers: any[]) {
  fn.updateTriggers = triggers;
  return fn;
}

export function withLoaderTriggers(fn: Accessor<any>, dataLoader: DataLoader) {
  fn.dataLoader = dataLoader;
  return withTriggers(fn, [dataLoader.id, dataLoader.updateTrigger]);
}
