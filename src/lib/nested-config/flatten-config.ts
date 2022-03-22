import { ConfigTree } from './config-tree';

export function flattenConfig<T>(configs: ConfigTree<T>) {
  const res: T[] = [];

  function flatten(configs: ConfigTree<T>) {
    for (const layer of configs) {
      if (Array.isArray(layer)) {
        flatten(layer);
      } else if (layer) {
        res.push(layer);
      }
    }
  }
  flatten(configs);

  return res;
}
