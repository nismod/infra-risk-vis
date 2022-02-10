import { DependencyList, EffectCallback, useEffect, useRef } from 'react';

// modified from https://stackoverflow.com/a/58125085/1478817
export const useChangeEffect = (
  /**
   * The callback to run when any of changeDeps changes.
   */
  effect: EffectCallback,

  /**
   * The full list of dependencies - equivalent to useEffect dependencies
   */
  deps: DependencyList,

  /**
   * The list of dependencies to watch for changes.
   */
  changeDeps: DependencyList,
) => {
  const changeRef = useRef(changeDeps || []);
  const initial = changeRef.current === changeDeps;
  const changeDepsUpdated = initial || !changeRef.current.every((w, i) => w === changeDeps[i]);
  changeRef.current = changeDeps;
  const nullDeps = deps.map(() => null);

  //eslint-disable-next-line react-hooks/exhaustive-deps
  return useEffect(changeDepsUpdated ? effect : () => {}, changeDepsUpdated ? deps : nullDeps);
};
