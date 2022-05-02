import { usePrevious } from 'lib/hooks/use-previous';
import { useEffect } from 'react';
import { RecoilState, useRecoilCallback, useRecoilValue } from 'recoil';
import { StateEffect } from './types';

/**
 * Run a state effect when a piece of state changes.
 * A state effect can modify other pieces of state.
 * @param state the recoil state to watch
 * @param effect the state effect to run when the state changes
 */
export function useStateEffect<T>(state: RecoilState<T>, effect: StateEffect<T>) {
  const stateValue = useRecoilValue(state);

  const previousStateValue = usePrevious(stateValue);

  const cb = useRecoilCallback(
    ({ transact_UNSTABLE }) =>
      (newValue: T, previousValue: T) => {
        transact_UNSTABLE((ops) => effect(ops, newValue, previousValue));
      },
    [effect],
  );

  useEffect(() => cb(stateValue, previousStateValue), [cb, stateValue, previousStateValue]);
}
