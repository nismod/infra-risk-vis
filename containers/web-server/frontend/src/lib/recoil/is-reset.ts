import { DefaultValue } from 'recoil';

export const isReset = (candidate: unknown): candidate is DefaultValue => {
  if (candidate instanceof DefaultValue) return true;
  return false;
};
