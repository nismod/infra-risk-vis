export type LayerTree<T> = (LayerTree<T> | false | null | undefined | T)[];

export function flattenLayers<T>(layers: LayerTree<T>) {
  const res: T[] = [];

  function flatten(layers: LayerTree<T>) {
    for (const layer of layers) {
      if (Array.isArray(layer)) {
        flatten(layer);
      } else if (layer) {
        res.push(layer);
      }
    }
  }
  flatten(layers);

  return res;
}
