import _ from 'lodash';
import { selectorFamily, waitForAll } from 'recoil';

import { RecoilReadableStateFamily } from './types';

export function groupedFamily<FVT, FPT>(
  key: string,
  family: RecoilReadableStateFamily<FVT, FPT>,
  paramsFamily: RecoilReadableStateFamily<string[], string>,
  paramFn: (group: string, param: string) => FPT,
) {
  return selectorFamily<Record<string, FVT>, string>({
    key,
    get:
      (group) =>
      ({ get }) => {
        const groupParams = get(paramsFamily(group));
        const deps = _.fromPairs(
          groupParams.map((param) => [param, family(paramFn(group, param))]),
        );
        return get(waitForAll(deps));
      },
  });
}
