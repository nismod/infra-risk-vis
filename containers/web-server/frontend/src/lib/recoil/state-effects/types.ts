import { TransactionInterface_UNSTABLE } from 'recoil';

export type StateEffectInterface = TransactionInterface_UNSTABLE;

/**
 * Type for a function to be used as a state effect.
 * The function will be called with the transaction interace,
 * the current value, and previous value of the tracked state.
 */
export type StateEffect<T> = (ops: StateEffectInterface, value: T, previousValue: T) => void;

/**
 * Like StateEffect<T> but without the previous value
 */
export type CurrentStateEffect<T> = (ops: StateEffectInterface, value: T) => void;
