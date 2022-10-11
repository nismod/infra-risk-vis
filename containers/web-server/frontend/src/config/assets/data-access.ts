import { dataLoaderManager } from 'lib/data-loader/data-loader-manager';
import { FieldSpec } from 'lib/data-map/view-layers';
import { extraProperty, featureProperty } from 'lib/deck/props/data-source';
import { withTriggers } from 'lib/deck/props/getters';
import { sumOrNone } from 'lib/helpers';

function getExpectedDamageKey(direct: boolean, hazard: string, rcp: string, epoch: number) {
  return `${direct ? 'ead' : 'eael'}__${hazard}__rcp_${rcp}__epoch_${epoch}__conf_None`;
}

const hazardTypes = ['fluvial', 'surface', 'coastal', 'cyclone', 'extreme_heat'];

function totalExpectedDamagesProperty(direct: boolean, { rcp, epoch }) {
  const hazardProperties = hazardTypes.map((ht) => featureProperty(getExpectedDamageKey(direct, ht, rcp, epoch)));

  return withTriggers((f) => sumOrNone(hazardProperties.map((p) => p(f))), [direct, rcp, epoch]);
}

export function getAssetDataAccessor(layer: string, fieldSpec: FieldSpec) {
  if (fieldSpec == null) return null;

  const { fieldGroup, fieldDimensions, field } = fieldSpec;

  if (fieldGroup === 'damages_expected') {
    const { hazard, rcp, epoch } = fieldDimensions;

    const isDirect = field.startsWith('ead_');

    if (hazard === 'all') {
      return totalExpectedDamagesProperty(isDirect, fieldDimensions);
    }
    return featureProperty(getExpectedDamageKey(isDirect, hazard, rcp, epoch));
  } else if (fieldGroup === 'damages_return_period') {
    // return return period damages dynamically loaded from API
    return extraProperty(dataLoaderManager.getDataLoader(layer, fieldSpec));
  } else if (fieldGroup === 'adaptation') {
    return extraProperty(dataLoaderManager.getDataLoader(layer, fieldSpec));
  } else {
    // field other than damages - use field name as key
    return featureProperty(field);
  }
}

export function assetDataAccessFunction(layer: string) {
  return (fieldSpec: FieldSpec) => getAssetDataAccessor(layer, fieldSpec);
}
