import _ from 'lodash';
import { atom, selector } from 'recoil';
import { recoilPersist } from 'recoil-persist';
import { NIL as NIL_UUID } from 'uuid';

import { Assessment } from 'config/assessment/assessment';
import { INTERVENTION_HIERARCHY } from 'config/assessment/interventions';
import { buildTreeConfig } from 'lib/controls/checkbox-tree/CheckboxTree';


const { persistAtom } = recoilPersist({
  key: 'infrastructure-resilience'
});

export const assessmentList = atom<Assessment[]>({
  key: 'assessmentList',
  default: [],
  effects: [persistAtom],
});

export const currentAssessmentID = atom<string>({
  key: 'currentAssessmentID',
  default: NIL_UUID,
});

export const currentAssessment = atom<Assessment>({
  key: 'currentAssessment',
  default: undefined
});

export const currentAssessmentInList = selector({
  key: 'currentAssessmentInList',
  get: ({get}) => {
    const assessmentID = get(currentAssessmentID);
    const assessments = get(assessmentList);
    return assessments.find((assessment) => assessment.id === assessmentID);
  },
  set: ({get, set}, newValue: Assessment) => {
    const assessments = get(assessmentList);
    const index = findIndexByID(assessments, newValue.id)
    console.log(index, newValue, assessments)
    index === -1?
      set(assessmentList, [...assessments, newValue])
      : set(assessmentList, replaceItemAtIndex(assessments, index, newValue))
  },
});

function findIndexByID(arr: Assessment[], id: string) {
  return arr.findIndex((item) => item.id === id);
}

function replaceItemAtIndex(arr: Assessment[], index: number, newValue: Assessment) {
  return [...arr.slice(0, index), newValue, ...arr.slice(index + 1)];
}

export const indicatorWeights = selector({
  key: 'indicatorWeights',
  get: ({ get }) => {
    const assessment = get(currentAssessment);
    return assessment.indicatorWeights;
  },
  set: ({ get, set }, newValue) => {
    const assessment = get(currentAssessment);
    set(currentAssessment, {
      ...assessment,
      indicatorWeights: newValue,
    });
  },
});

export const interventionSelection = selector({
  key: 'interventionSelection',
  get: ({ get }) => {
    const assessment = get(currentAssessment);
    return assessment.interventionSelection;
  },
  set: ({ get, set }, newValue) => {
    const assessment = get(currentAssessment);
    set(currentAssessment, {
      ...assessment,
      interventionSelection: newValue,
    });
  },
});

export const interventionTreeConfig = buildTreeConfig(INTERVENTION_HIERARCHY);