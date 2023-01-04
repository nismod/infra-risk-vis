import _ from 'lodash';

export type ParamValue = any;
export type ParamDomain<PT extends ParamValue = ParamValue> = PT[];

export type ParamGroup = Record<string, ParamValue>;

export type ParamDependency<PT extends ParamValue, PGT extends ParamGroup> = (
  params: PGT,
) => ParamDomain<PT>;

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
    const newParamOptions = getNewParamOptions(
      resolvedParams,
      paramDomains[param],
      paramDependencies[param],
    );

    // if the new options don't include the current param value, switch value to the first option
    if (!newParamOptions.includes(paramValue)) {
      resolvedParams[param as K] = newParamOptions[0];
    }

    newOptions[param as K] = newParamOptions;
  }
  return [resolvedParams, newOptions];
}

/*
 * Inferring domains / dependencies from data
 */

type DependenciesSpec<T extends object> = Record<keyof T, (keyof T)[]>;

type SortIteratees<T> = _.Many<_.ListIteratee<T>>; // taken from lodash _.sortBy
type SortSpec<T extends object> = Partial<Record<keyof T, SortIteratees<any>>>;

function getGroupKey<T>(keys: (keyof T)[]) {
  return (obj: T) => keys.map((k) => obj[k]).join('+');
}

/**
 * group a collection by multiple properties
 * @param data the collection to group
 * @param properties list of string names of the properties
 * @returns dictionary of arrays, keyed by a generated group key
 */
function groupByMulti<T>(data: T[], properties: (keyof T)[]) {
  return _.groupBy(data, getGroupKey(properties));
}

function makeDependencyFunction<T extends object>(
  data: T[],
  param: keyof T,
  dependencies: (keyof T)[],
  sortByIteratees: SortIteratees<T[typeof param]> = _.identity,
) {
  // group the data objected by co-occurring values of the dependency fields
  const grouped = groupByMulti(data, dependencies);

  // for each combination of dependency values, get a sorted and unique list
  // of values of the dependent parameter
  const groupedDomains = _.mapValues(grouped, (objects) =>
    _(objects)
      .map((obj) => obj[param])
      .uniq()
      .sortBy(sortByIteratees)
      .value(),
  );

  // return a function that gets the possible options of the dependent parameter based on a current state of the params
  return (params: T) => groupedDomains[getGroupKey(dependencies)(params)];
}

export function inferDependenciesFromData<T extends object>(
  data: T[],
  depSpec: DependenciesSpec<T>,
  sortSpec?: SortSpec<T>,
): ParamGroupDependencies<T> {
  const dependencies: ParamGroupDependencies<T> = {};
  for (const [param, inputs] of Object.entries(depSpec) as [keyof T, (keyof T)[]][]) {
    if (inputs.length > 0) {
      dependencies[param] = makeDependencyFunction(
        data,
        param as keyof T,
        inputs as (keyof T)[],
        sortSpec?.[param],
      );
    }
  }
  return dependencies;
}

export function inferDomainsFromData<T extends object>(data: T[]): ParamGroupDomains<T> {
  const domains: any = {};
  const keys = Object.keys(data[0]);
  for (const key of keys) {
    const values = data.map((d) => d[key]);
    domains[key] = _.uniq(values);
  }
  return domains;
}
