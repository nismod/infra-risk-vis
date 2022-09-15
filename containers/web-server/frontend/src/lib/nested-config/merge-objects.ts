type MergeStrategy = <T>(oldValue: T, newValue: T) => T;

/**
 * A function to merge multiple props-like objects.
 * This is inspired by Deck.GL layer props-merging behavior, but is extended
 * with the ability to specify custom merging strategies for entries specified by key
 * @param objects
 * @param mergeStrategies
 * @returns
 */
export function mergeObjects(objects: object[], mergeStrategies: Record<string, MergeStrategy>): object {
  const mergedProps = {};

  for (const props of objects) {
    for (const [key, value] of Object.entries(props)) {
      mergedProps[key] = mergeStrategies[key] ? mergeStrategies[key](mergedProps[key], value) : value;
    }
  }
  return mergedProps;
}

export const mergeValue = (oldVal, newVal) => Object.assign({}, oldVal ?? {}, newVal);
