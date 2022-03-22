import { MVTLayer, TileLayer, BitmapLayer, GeoJsonLayer } from 'deck.gl';

import { mergeDeckProps } from './merge-props';

export function mvtLayer(...props) {
  return new MVTLayer(mergeDeckProps(props));
}

export function tileLayer(...props) {
  return new TileLayer(mergeDeckProps(props));
}

export function bitmapLayer(...props) {
  return new BitmapLayer(mergeDeckProps(props));
}

export function geoJsonLayer(...props) {
  return new GeoJsonLayer(mergeDeckProps(props));
}
