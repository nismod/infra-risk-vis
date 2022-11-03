import { DataFilterExtension } from '@deck.gl/extensions/typed';

/**
 * Filter features by ID (or unique ID property) using the GPU extension
 * @param featureId the feature ID to filter by
 * @param uniqueIdProperty the name of a unique ID property, if feature IDs are not present
 * @returns deck.gl layer props that configure feature filtering
 */
export function featureFilter(featureId: string | number, uniqueIdProperty?: string) {
  const filterFn = uniqueIdProperty
    ? (x) => (x.properties[uniqueIdProperty] === featureId ? 1 : 0)
    : (x) => (x.id === featureId ? 1 : 0);

  return {
    updateTriggers: {
      getFilterValue: [uniqueIdProperty, featureId],
    },

    getFilterValue: filterFn,
    filterRange: [1, 1],
    extensions: [new DataFilterExtension({ filterSize: 1 })],
  };
}
