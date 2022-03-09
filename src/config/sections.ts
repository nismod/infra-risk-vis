import { StyleSelectionOption } from 'state/sections';
import { BUILDING_STYLES } from './buildings/styles';
import { NETWORK_STYLES } from './networks/styles';
import { REGION_STYLES } from './regions/styles';

export const SECTIONS_CONFIG: Record<string, { styles?: Record<string, StyleSelectionOption> }> = {
  assets: {
    styles: NETWORK_STYLES,
  },
  hazards: {},
  buildings: {
    styles: BUILDING_STYLES,
  },
  regions: {
    styles: REGION_STYLES,
  },
};