export type Param = any;
export type ParamDomain<PT extends Param = Param> = PT[];

export type ParamGroup = Record<string, Param>;

export type ParamDependency<PT extends Param, PGT extends ParamGroup> = (params: PGT) => ParamDomain<PT>;

export type ParamGroupDomains<PGT extends ParamGroup = ParamGroup> = {
  [K in keyof PGT]: ParamDomain<PGT[K]>;
};

export type ParamGroupDependencies<PGT extends ParamGroup = ParamGroup> = {
  [K in keyof PGT]?: ParamDependency<PGT[K], PGT>;
};

export interface DataParamGroupConfig<PGT extends ParamGroup = ParamGroup> {
  paramDomains: ParamGroupDomains<PGT>;
  paramDefaults: PGT;
  paramDependencies?: ParamGroupDependencies<PGT>;
}

function getNewParamOptions(updatedParams, domain, dependencyFn) {
  return dependencyFn?.(updatedParams) ?? domain;
}

export function resolveParamDependencies<PGT extends ParamGroup = ParamGroup>(
  updatedParams: PGT,
  config: DataParamGroupConfig<PGT>,
): [PGT, ParamGroupDomains<PGT>] {
  type K = keyof PGT;
  const { paramDomains, paramDependencies = {} } = config;

  const resolvedParams = { ...updatedParams };
  const newOptions: ParamGroupDomains<PGT> = {} as any;

  for (const [param, paramValue] of Object.entries(updatedParams)) {
    const newParamOptions = getNewParamOptions(updatedParams, paramDomains[param], paramDependencies[param]);

    // if the new options don't include the current param value, switch value to the first option
    if (!newParamOptions.includes(paramValue)) {
      resolvedParams[param as K] = newParamOptions[0];
    }

    newOptions[param as K] = newParamOptions;
  }
  return [resolvedParams, newOptions];
}
