import { ValueLabel } from '@/lib/controls/params/value-label';
import { makeColorConfig } from '@/lib/helpers';

export const PROTECTED_AREA_TYPES = ['land', 'marine'] as const;

export type ProtectedAreaType = typeof PROTECTED_AREA_TYPES[number];

export const PROTECTED_AREA_LABELS: ValueLabel<ProtectedAreaType>[] = [
  {
    value: 'land',
    label: 'Terrestrial and Inland Waters',
  },
  {
    value: 'marine',
    label: 'Marine',
  },
];

export const PROTECTED_AREA_COLORS = makeColorConfig<ProtectedAreaType>({
  marine: '#004DA8',
  land: '#38A800',
});
