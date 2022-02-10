import { RecoilState } from 'recoil';

export type RecoilStateFamily<DataType, ParamType> = (param: ParamType) => RecoilState<DataType>;
