import { TransactionInterface_UNSTABLE } from 'recoil';

export type StateEffect<T> = (ops: TransactionInterface_UNSTABLE, value: T) => void;
