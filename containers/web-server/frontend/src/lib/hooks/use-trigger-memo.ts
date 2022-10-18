import { DependencyList, useMemo } from 'react';

/**
 * Standard useMemo that accepts an additional trigger that is not a dependency of the callback,
 * and which, when updated, will cause a recalculation.
 * This is mostly to avoid having to manually silence rules-of-hooks warnings about unnecessary deps,
 * when the trigger is not used in the callback
 * @param callback function to call to return the memoed result
 * @param dependencies dependencies of the callback
 * @param trigger trigger which will cause recalculation. It shouldn't be a dependency of the callback
 * @returns memoed result of the callback
 */
export function useTriggerMemo<T>(callback: () => T, dependencies: DependencyList, trigger: any) {
  return useMemo(
    callback,
    // eslint-disable-next-line react-hooks/exhaustive-deps
    [callback, ...dependencies, trigger],
  );
}
