import _ from 'lodash';
import { RecoilValue, atomFamily } from 'recoil';

export const hazardSelectionState = atomFamily({
  key: 'hazardSelectionState',
  default: false,
});

interface TransactionGetterInterface {
  get<T>(a: RecoilValue<T>): T;
}

export function getHazardSelectionAggregate({ get }: TransactionGetterInterface, hazards: string[]) {
  return _.fromPairs(hazards.map((group) => [group, get(hazardSelectionState(group))]));
}
