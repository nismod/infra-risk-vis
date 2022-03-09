import { SECTIONS_CONFIG } from 'config/sections';
import { VIEW_SECTIONS } from 'config/views';
import { isReset } from 'lib/recoil/is-reset';
import _ from 'lodash';
import { atom, selector } from 'recoil';

import {
  sectionStyleOptionsState,
  sectionStyleValueState,
  sectionVisibilityState,
  sidebarSectionExpandedState,
} from 'state/sections';

export const viewImpl = atom({
  key: 'viewImpl',
  default: 'exposure',
});

export const viewState = selector<string>({
  key: 'viewState',
  get: ({ get }) => get(viewImpl),
  set: ({ set }, view) => {
    // ignore reset of view - potentially unsafe?
    if (isReset(view)) return;

    set(viewImpl, view);

    const viewSectionsConfig = VIEW_SECTIONS[view];

    _.forEach(viewSectionsConfig, (sectionConfig, section) => {
      set(sectionVisibilityState(section), sectionConfig.visible);
      set(sidebarSectionExpandedState(section), sectionConfig.expanded);
      const styleOptions = sectionConfig.styles?.map((style) => SECTIONS_CONFIG[section].styles[style]);
      set(sectionStyleOptionsState(section), styleOptions);
      set(sectionStyleValueState(section), sectionConfig.defaultStyle);
    });
  },
});
