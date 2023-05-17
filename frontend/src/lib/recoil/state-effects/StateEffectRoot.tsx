import { FC } from 'react';
import { useStateEffect } from './use-state-effect';

/**
 * A component that wraps a state effect to prevent a parent from being unnecessarily re-rendered
 */
export const StateEffectRoot: FC<{ state; effect }> = ({ state, effect }) => {
  useStateEffect(state, effect);
  return null;
};
