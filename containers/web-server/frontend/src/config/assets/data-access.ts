import { dataLoaderManager } from '@/lib/data-loader/data-loader-manager';
import { FieldSpec } from '@/lib/data-map/view-layers';
import { extraProperty, featureProperty } from '@/lib/deck/props/data-source';
import { withTriggers } from '@/lib/deck/props/getters';
import { sumOrNone } from '@/lib/helpers';

import { HAZARD_TYPES, HazardType } from '../hazards/metadata';

/**
 * handle quirks/inconsistencies in feature property naming for hazards
 */
function lookupHazard(hazard: string) {
  if (hazard === 'fluvial') return 'river';
  return hazard;
}

/**
 * handle quirks/inconsistencies in feature property naming for RCP
 */
function lookupRcp(rcp: string, layer: string) {
  if (rcp == null) return rcp;
  if (rcp === 'baseline') {
    if (layer.startsWith('road_')) {
      return 'historical';
    }
    return rcp;
  }
  return `rcp${rcp.replace(/\./g, 'p')}`;
}

/**
 * handle quirks/inconsistencies in feature property naming for epoch
 */
function lookupEpoch(epoch: string, hazard: HazardType, layer: string) {
  if (epoch === 'present') {
    if (hazard === 'cyclone') {
      return '2020';
    }
    if (layer.startsWith('road_')) {
      return '1980';
    }
  }
  return epoch;
}

/**
 * Construct the feature property name based on what data is requested
 * Handle various quirks in the data naming
 */
function getExpectedDamageKey(layer: string, direct: boolean, hazard: HazardType, rcp: string, epoch: string) {
  return `${direct ? 'ead' : 'eael'}__${lookupHazard(hazard)}__rcp_${lookupRcp(rcp, layer)}__epoch_${lookupEpoch(
    epoch,
    hazard,
    layer,
  )}`;
}

/**
 * Make a value accessor that sums up individual hazard damages on the fly into a total figure
 * NOTE: this isn't currently used in the G-SRAT version of the tool
 * The code here assumes damages for all hazards contribute to a total figure
 */
function totalExpectedDamagesProperty(layer: string, direct: boolean, { rcp, epoch }) {
  const hazardProperties = HAZARD_TYPES.map((ht) =>
    featureProperty(getExpectedDamageKey(layer, direct, ht, rcp, epoch)),
  );

  return withTriggers((f) => sumOrNone(hazardProperties.map((p) => p(f))), [direct, rcp, epoch]);
}

/**
 * Defines how to get data for an asset layer, depending on what field is required
 */
export function getAssetDataAccessor(layer: string, fieldSpec: FieldSpec) {
  if (fieldSpec == null) return null;

  const { fieldGroup, fieldDimensions, field } = fieldSpec;

  if (fieldGroup === 'damages_expected') {
    /**
     * expected damages are stored in vector feature properties - need
     */
    const { hazard, rcp, epoch } = fieldDimensions;

    const isDirect = field.startsWith('ead_');

    /**
     * total hazard damages need to be calculated on the fly from the individual properties
     */
    if (hazard === 'all') {
      return totalExpectedDamagesProperty(layer, isDirect, fieldDimensions);
    }
    return featureProperty(getExpectedDamageKey(layer, isDirect, hazard, rcp, epoch));
  } else if (fieldGroup === 'damages_return_period') {
    /**
     * return period damages are loaded dynamically from the features API
     */
    return extraProperty(dataLoaderManager.getDataLoader(layer, fieldSpec));
  } else if (fieldGroup === 'adaptation') {
    /**
     * adaptation option values are loaded dynamically from the features API
     */
    return extraProperty(dataLoaderManager.getDataLoader(layer, fieldSpec));
  } else {
    /**
     * for any other fields, assume it's a property contained in the vector features
     */
    return featureProperty(field);
  }
}

export function assetDataAccessFunction(layer: string) {
  return (fieldSpec: FieldSpec) => getAssetDataAccessor(layer, fieldSpec);
}
