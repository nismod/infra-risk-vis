import { RecoilValueReadOnly, atom, selector } from 'recoil';

function firstElem<T>(options: T) {
  return options?.[0];
}

export type DefaultFunction<T> = (options: T[]) => T;

export function makeSelectState<T>(
  key: string,
  optionsState: RecoilValueReadOnly<T[]>,
  allowEmpty: boolean = false,
  defaultFn: DefaultFunction<T> = firstElem,
) {
  const selectedImpl = atom({
    key: `${key}/impl`,
    default: null,
  });

  const defaultImpl = selector({
    key: `${key}/default`,
    get: ({ get }) => defaultFn(get(optionsState)),
  });

  const resultState = selector({
    key,
    get: ({ get }) => {
      const selectedOption = get(selectedImpl);
      const options = get(optionsState);
      if ((selectedOption == null && allowEmpty) || !options.includes(selectedOption)) {
        return get(defaultImpl);
      }
      return selectedOption;
    },
    set: ({ set }, newValue) => set(selectedImpl, newValue),
  });

  return resultState;
}
