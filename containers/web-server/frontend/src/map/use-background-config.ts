import _ from 'lodash';
import { useMemo } from 'react';

import { BACKGROUNDS, BackgroundName } from '@/config/backgrounds';

function visible(isVisible: boolean): 'visible' | 'none' {
  return isVisible ? 'visible' : 'none';
}

function makeBackgroundConfig(background: BackgroundName) {
  return {
    version: 8,
    sources: _.mapValues(BACKGROUNDS, (b) => b.source),
    layers: Object.values(BACKGROUNDS).map((b) =>
      _.merge(b.layer, { layout: { visibility: visible(background === b.id) } }),
    ),
  };
}

export function useBackgroundConfig(background: BackgroundName) {
  return useMemo(() => makeBackgroundConfig(background), [background]);
}
