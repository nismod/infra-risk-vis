import { atomFamily } from 'recoil';

export const sectionVisibilityState = atomFamily<boolean, string>({
  key: 'sectionVisibilityState',
  default: true,
});

export const sidebarSectionExpandedState = atomFamily({
  key: 'sidebarSectionExpandedState',
  default: false,
});

export const sectionStyleValueState = atomFamily<string, string>({
  key: 'sectionStyleValueState',
  default: '',
});

export interface StyleSelectionOption {
  id: string;
  label: string;
}

export const sectionStyleOptionsState = atomFamily<StyleSelectionOption[], string>({
  key: 'sectionStyleOptionsState',
  default: [],
});

export const sectionStyleDefaultValueState = atomFamily<string, string>({
  key: 'sectionStyleDefaultValueState',
  default: null,
});
