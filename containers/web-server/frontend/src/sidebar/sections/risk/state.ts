import { atom } from 'recoil';

import { ValueLabel } from '@/lib/controls/params/value-label';
import { StateEffect } from '@/lib/recoil/state-effects/types';

import { sidebarVisibilityToggleState } from '@/sidebar/SidebarContent';

export const RISK_TARGETS = ['population', 'buildings', 'transport', 'power'] as const;
export type RiskTargetType = typeof RISK_TARGETS[number];

export const RISK_TARGET_OPTIONS: ValueLabel<RiskTargetType>[] = [
  { value: 'population', label: 'Population' },
  { value: 'buildings', label: 'Buildings' },
  { value: 'transport', label: 'Roads' },
  { value: 'power', label: 'Power Grid' },
];

export const riskTargetTypeState = atom<RiskTargetType>({
  key: 'riskTargetTypeState',
  default: 'population',
});

const RISK_TARGET_SIDEBAR_SECTIONS: Record<RiskTargetType, string> = {
  population: 'exposure/population',
  buildings: 'exposure/buildings',
  power: 'exposure/infrastructure',
  transport: 'exposure/infrastructure',
};

export const riskTargetTypeStateEffect: StateEffect<RiskTargetType> = ({ set }, riskTarget, previousRiskTarget) => {
  const previousTargetSection = RISK_TARGET_SIDEBAR_SECTIONS[previousRiskTarget];
  const newTargetSection = RISK_TARGET_SIDEBAR_SECTIONS[riskTarget];

  set(sidebarVisibilityToggleState(previousTargetSection), false);
  set(sidebarVisibilityToggleState(newTargetSection), true);
};

export const RISK_SOURCES = ['fluvial', 'cyclone', 'extreme_heat', 'drought', 'earthquake'] as const;

export type RiskSourceType = typeof RISK_SOURCES[number];

export const RISK_SOURCE_OPTIONS: ValueLabel<RiskSourceType>[] = [
  { value: 'fluvial', label: 'River Flooding' },
  { value: 'cyclone', label: 'Tropical Cyclones' },
  { value: 'extreme_heat', label: 'Extreme Heat' },
  { value: 'drought', label: 'Droughts' },
  { value: 'earthquake', label: 'Seismic Risk' },
];

export const riskSourceTypeState = atom<RiskSourceType>({
  key: 'riskSourceTypeState',
  default: 'fluvial',
});
