import { dataLoaderManager } from 'lib/data-loader/data-loader-manager';
import { FieldSpec } from 'lib/data-map/view-layers';
import { extraProperty, featureProperty } from 'lib/deck/props/data-source';
import { withTriggers } from 'lib/deck/props/getters';
import { sumOrNone } from 'lib/helpers';

function getDamageKey(hazard: string, rcp: string, epoch: number) {
  return `${hazard}__rcp_${rcp}__epoch_${epoch}__conf_None`;
}

const hazardTypes = ['fluvial', 'surface', 'coastal', 'cyclone'];

function totalDamagesProperty({ rcp, epoch }) {
  const hazardProperties = hazardTypes.map((ht) => featureProperty(getDamageKey(ht, rcp, epoch)));

  return withTriggers((f) => sumOrNone(hazardProperties.map((p) => p(f))), [rcp, epoch]);
}

export function getAssetDataAccessor(layer: string, fieldSpec: FieldSpec) {
  const { field, fieldParams } = fieldSpec;

  if (field === 'damages') {
    const { damage_type, hazard, rcp, epoch } = fieldParams;

    if (damage_type === 'indirect') {
      // load indirect damages from API
      return extraProperty(dataLoaderManager.getDataLoader(layer, fieldSpec));
    }
    if (hazard === 'total-damages') {
      return totalDamagesProperty(fieldParams);
    }
    return featureProperty(getDamageKey(hazard, rcp, epoch));
  } else {
    // field other than damages - use field name as key
    return featureProperty(field);
  }
}

export function assetDataAccessFunction(layer: string) {
  return ({ styleParams }) => {
    if (styleParams?.colorMap) {
      const { colorField } = styleParams.colorMap;

      return {
        dataAccessor: getAssetDataAccessor(layer, colorField),
        dataLoader: dataLoaderManager.getDataLoader(layer, colorField),
      };
    }
  };
}
