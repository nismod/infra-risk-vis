import { atom, atomFamily } from 'recoil';
import { InteractionTarget } from './use-interactions';

export function hasHover(target: InteractionTarget<any> | InteractionTarget<any>[]) {
  if (Array.isArray(target)) {
    return target.length > 0;
  }
  return !!target;
}

export const hoverState = atomFamily<InteractionTarget<any> | InteractionTarget<any>[], string>({
  key: 'hoverState',
  default: null,
});

export const hoverPositionState = atom({
  key: 'hoverPosition',
  default: null,
});

export const selectionState = atomFamily<InteractionTarget<any>, string>({
  key: 'selectionState',
  default: null,
});
