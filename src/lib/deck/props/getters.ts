export type Accessor<To, From = any> = ((x: From) => To) & { updateTriggers?: any[] };
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
