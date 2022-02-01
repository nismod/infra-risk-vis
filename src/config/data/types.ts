type Params = object;
type ParamDomain<PT> = PT[];

export type ParamDomains<PT extends Params> = {
  [K in keyof PT]: ParamDomain<PT[K]>;
};

export type ParamDependencies<PT extends Params> = {
  [K in keyof PT]?: (params: PT) => ParamDomain<PT[K]>;
};

export interface DataParamConfig<P extends Params> {
  paramDomains: ParamDomains<P>;
  paramDefaults: P;
  paramDependencies?: ParamDependencies<P>;
}
