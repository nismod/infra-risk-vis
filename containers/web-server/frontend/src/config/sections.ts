import { StyleSelectionOption } from '@/state/sections';

import { NETWORK_STYLES } from './networks/styles';

export const SECTIONS_CONFIG: Record<string, { styles?: Record<string, StyleSelectionOption> }> = {
  assets: {
    styles: NETWORK_STYLES,
  },
  hazards: {},
  buildings: {},
};
