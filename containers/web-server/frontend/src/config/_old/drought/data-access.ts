import { FieldSpec } from '@/lib/data-map/view-layers';
import { featureProperty } from '@/lib/deck/props/data-source';

function getDroughtPropertyKey(field: string, rcp?: string) {
  return `${field}${rcp ? `__rcp_${rcp}` : ''}`;
}

export function getDroughtDataAccessor(fieldSpec: FieldSpec) {
  if (fieldSpec == null) return null;

  const { field, fieldDimensions } = fieldSpec;

  return featureProperty(getDroughtPropertyKey(field, fieldDimensions.rcp));
}
