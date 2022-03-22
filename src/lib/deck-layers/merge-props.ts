import { ConfigTree, flattenConfig } from 'lib/config-tree';

/**
 * A function to merge multiple props objects passed to a Deck.GL layer.
 * This extends the base Deck.GL behaviour in a few ways:
 * - falsy elements of the array are ignored
 * - nested arrays are flattened
 * - compound props (currently only `updateTriggers`) are merged instead of overwritten
 * @param propsArray
 */
export const mergeDeckProps = (...propsArray: ConfigTree<object>) => {
  const flattenedProps = flattenConfig(propsArray);
  const mergedProps: Record<string, any> = {};

  for (const props of flattenedProps) {
    for (const [key, value] of Object.entries(props)) {
      let newVal: any;
      if (key === 'updateTriggers') {
        const oldVal = mergedProps[key] ?? {};
        newVal = Object.assign({}, oldVal, value);
      } else {
        newVal = value;
      }
      mergedProps[key] = newVal;
    }
  }
  return mergedProps;
};
