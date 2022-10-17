import { BitmapLayer, GeoJsonLayer, MVTLayer, TileLayer } from 'deck.gl';

import { ConfigTree } from '@/lib/nested-config/config-tree';
import { flattenConfig } from '@/lib/nested-config/flatten-config';
import { mergeObjects, mergeValue } from '@/lib/nested-config/merge-objects';

const deckPropsMergeStrategies = {
  updateTriggers: mergeValue,
};

/**
 * A function to merge multiple props objects passed to a Deck.GL layer.
 * This extends the base Deck.GL behaviour in a few ways:
 * - falsy elements of the array are ignored
 * - nested arrays are flattened
 * - compound props (currently only `updateTriggers`) are merged instead of overwritten
 */
function mergeDeckProps(...props: ConfigTree<object>): any {
  const flattenedProps = flattenConfig(props);

  return mergeObjects(flattenedProps, deckPropsMergeStrategies);
}

function wrap(deckClass) {
  return (...props) => new deckClass(mergeDeckProps(props));
}

export const mvtLayer = wrap(MVTLayer);
export const tileLayer = wrap(TileLayer);
export const bitmapLayer = wrap(BitmapLayer);
export const geoJsonLayer = wrap(GeoJsonLayer);
