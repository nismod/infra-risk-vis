import { dataLoaderManager } from 'lib/data-loader/data-loader-manager';
import { FieldSpec } from 'lib/data-map/view-layers';
import { extraProperty, featureProperty } from 'lib/deck/props/data-source';
import { withTriggers } from 'lib/deck/props/getters';
import { sumOrNone } from 'lib/helpers';

function getExpectedDirectDamageKey(isDirect: boolean, hazard: string, rcp: string, epoch: number) {
  return `${isDirect ? 'ead' : 'eael'}__${hazard}__rcp_${rcp}__epoch_${epoch}__conf_None`;
}

const hazardTypes = ['fluvial', 'surface', 'coastal', 'cyclone'];

function totalExpectedDirectDamagesProperty(isDirect: boolean, { rcp, epoch }) {
  const hazardProperties = hazardTypes.map((ht) =>
    featureProperty(getExpectedDirectDamageKey(isDirect, ht, rcp, epoch)),
  );

  return withTriggers((f) => sumOrNone(hazardProperties.map((p) => p(f))), [rcp, epoch]);
}

export function getAssetDataAccessor(layer: string, fieldSpec: FieldSpec) {
  const { fieldGroup, fieldDimensions, field } = fieldSpec;

  if (fieldGroup === 'damages_expected') {
    const { hazard, rcp, epoch } = fieldDimensions;

    const isDirect = field.startsWith('ead_');

    if (hazard === 'all') {
      return totalExpectedDirectDamagesProperty(isDirect, fieldDimensions);
    }
    return featureProperty(getExpectedDirectDamageKey(isDirect, hazard, rcp, epoch));
  } else if (fieldGroup === 'damages_return_period') {
    // return return period damages dynamically loaded from API
    return extraProperty(dataLoaderManager.getDataLoader(layer, fieldSpec));
  } else {
    // field other than damages - use field name as key
    return featureProperty(field);
  }
}

export function assetDataAccessFunction(layer: string) {
  return ({ styleParams }) => {
    if (styleParams?.colorMap) {
      const { colorField } = styleParams.colorMap;

      const accessor = getAssetDataAccessor(layer, colorField);
      return {
        dataAccessor: accessor,
        dataLoader: accessor.dataLoader,
      };
    }
  };
}
