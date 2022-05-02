import { StyleSelectionOption } from 'state/sections';
import { BUILDING_STYLES } from './buildings/styles';
import { DROUGHT_STYLES } from './drought/styles';
import { NETWORK_STYLES } from './networks/styles';
import { REGION_STYLES } from './regions/styles';
import { MARINE_STYLES, TERRESTRIAL_STYLES } from './solutions/styles';

export const SECTIONS_CONFIG: Record<string, { styles?: Record<string, StyleSelectionOption> }> = {
  assets: {
    styles: NETWORK_STYLES,
  },
  drought: {
    styles: DROUGHT_STYLES,
  },
  hazards: {},
  buildings: {
    styles: BUILDING_STYLES,
  },
  regions: {
    styles: REGION_STYLES,
  },
  terrestrial: {
    styles: TERRESTRIAL_STYLES,
  },
  marine: {
    styles: MARINE_STYLES,
  },
};
