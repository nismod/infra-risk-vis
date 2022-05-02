export interface ValueLabel<K extends string | number = string> {
  value: K;
  label: string;
}

export function isValueLabel(value): value is ValueLabel {
  return typeof value === 'object' && 'value' in value && 'label' in value;
}

export function getValueLabel(value: any): ValueLabel {
  if (isValueLabel(value)) {
    return value;
  }

  return {
    value: value,
    label: value,
  };
}
