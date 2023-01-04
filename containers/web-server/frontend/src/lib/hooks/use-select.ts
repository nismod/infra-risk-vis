import { useMemo, useState } from 'react';

export type DefaultFunction<T> = (options: T[]) => T;

function firstElem<T>(options: T) {
  return options?.[0];
}

export function useSelect<T>(
  options: T[],
  allowEmpty: boolean = false,
  defaultFn: DefaultFunction<T> = firstElem,
) {
  const [selectedOption, setSelectedOption] = useState(defaultFn(options));

  const defaultValue = useMemo(() => defaultFn(options), [options, defaultFn]);

  const resultSelection = useMemo(() => {
    if ((selectedOption == null && allowEmpty) || !options.includes(selectedOption)) {
      return null;
    }
    return selectedOption;
  }, [selectedOption, options, allowEmpty]);

  return [resultSelection ?? defaultValue, setSelectedOption] as const;
}
