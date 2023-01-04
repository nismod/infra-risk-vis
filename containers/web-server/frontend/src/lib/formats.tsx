import { ReactNode } from 'react';

export type FormatFunction<T = any> = (x: T) => ReactNode | string;

/**
 * Helper function to format numerical values in legends, tooltips etc.
 * TODO: make usage of this function and numFormat etc consistent
 * @param format simple format string where _ will be replaced by the formatted number
 * @param numberFormatOptions options to pass to Number.toLocaleString()
 * @returns function that accepts a number and returns a formatted string
 */
export function makeValueFormat(
  format: string | ((children: string) => ReactNode),
  numberFormatOptions?: Intl.NumberFormatOptions,
): FormatFunction<number> {
  const formatFn =
    typeof format === 'string' ? (x) => <>{format.replace('_', x)}</> : (x) => format(x);
  return (value: number) => formatFn(value.toLocaleString(undefined, numberFormatOptions));
}

/**
 * Wraps a format function in a null/undefined check and returns
 * a replacement string ('-' by default) if the value is null
 */
export function nullFormat<T>(
  formatFn: FormatFunction<T>,
  nullReplacement: string = '-',
): FormatFunction<T> {
  return (x: any) => (x != null ? formatFn(x) : nullReplacement);
}
