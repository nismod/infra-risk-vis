import { RecoilState, RecoilValueReadOnly } from 'recoil';

/**
 * A readable and writable recoil state family
 */
export type RecoilStateFamily<DataType, ParamType> = (param: ParamType) => RecoilState<DataType>;

/**
 * A recoil state family for which only read operation is required
 */
export type RecoilReadableStateFamily<DataType, ParamType> = (
  param: ParamType,
) => RecoilState<DataType> | RecoilValueReadOnly<DataType>;
