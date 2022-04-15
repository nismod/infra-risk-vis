import * as d3 from 'd3-color';
import _ from 'lodash';

/**
 * Common helper functions
 *
 */

export function commas(x) {
  return x.toString().replace(/\B(?=(\d{3})+(?!\d))/g, ',');
}

export function titleCase(str: string) {
  var splitStr = str.toLowerCase().split(' ');
  for (var k = 0; k < splitStr.length; k++) {
    splitStr[k] = splitStr[k].charAt(0).toUpperCase() + splitStr[k].substring(1);
  }
  return splitStr.join(' ');
}

const FORMATTER = Intl.NumberFormat('en-GB', { maximumSignificantDigits: 3 });

export function numFormat(n: number) {
  return FORMATTER.format(n);
}

export function unique<T>(arr: T[]) {
  return Array.from(new Set(arr));
}

export function colorCssToRgb(cssColor: string): [number, number, number, number?] {
  const { r, g, b } = d3.rgb(cssColor);
  return [r, g, b];
}

export function toDictionary<T, K extends string, V>(
  array: T[],
  keyFn: (x: T) => K,
  valueFn: (x: T) => V,
): Record<K, V> {
  return Object.fromEntries(array.map((x) => [keyFn(x), valueFn(x)])) as Record<K, V>;
}

export function makeConfig<C, K extends string>(cfg: (C & { id: K })[]) {
  return toDictionary(
    cfg,
    (x) => x.id,
    (x) => x,
  );
}

export function makeColorConfig<K extends string>(cfg: Record<K, string>) {
  return _.mapValues(cfg, (c) => ({ css: c, deck: colorCssToRgb(c) }));
}

// see discussion at https://stackoverflow.com/questions/23437476/in-typescript-how-to-check-if-a-string-is-numeric
export function isNumeric(val: any): boolean {
  return !(val instanceof Array) && val - parseFloat(val) + 1 >= 0;
}

export function truthyKeys(obj: Record<string, any>) {
  return Object.keys(obj).filter((k) => obj[k]);
}

/**
 * Returns sum of array elements, or null if all elements are null
 */
export function sumOrNone(arr: number[]): number | null {
  let result: number = null;

  for (const x of arr) {
    if (x != null) {
      if (result == null) result = x;
      else result += x;
    }
  }
  return result;
}

// adapted from https://stackoverflow.com/a/44661948/1478817
export function downloadFile(content: string, mimeType: string, fileName: string) {
  const element = document.createElement('a');
  const file = new Blob([content], { type: mimeType });
  element.href = URL.createObjectURL(file);
  element.download = fileName;
  document.body.appendChild(element); // Required for this to work in FireFox
  element.click();
}
