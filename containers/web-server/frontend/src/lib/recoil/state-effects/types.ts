import { TransactionInterface_UNSTABLE } from 'recoil';

export type StateEffectInterface = TransactionInterface_UNSTABLE;
export type StateEffect<T> = (ops: StateEffectInterface, value: T, previousValue: T) => void;
