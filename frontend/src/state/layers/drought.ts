import { ASSETS_SOURCE } from 'config/assets/source';
import { VECTOR_COLOR_MAPS } from 'config/color-maps';
import { getDroughtDataAccessor } from 'config/drought/data-access';
import { getDroughtOptionsDataFormats, getDroughtRiskDataFormats } from 'config/drought/data-formats';
import {
  DroughtOptionsVariableType,
  DroughtRiskVariableType,
  DROUGHT_OPTIONS_VARIABLES_WITH_RCP,
  DROUGHT_RISK_VARIABLES_WITH_RCP,
} from 'config/drought/metadata';
import { colorMap } from 'lib/color-map';
import { ColorSpec, FieldSpec, ViewLayer } from 'lib/data-map/view-layers';
import { selectableMvtLayer } from 'lib/deck/layers/selectable-mvt-layer';
import { dataColorMap } from 'lib/deck/props/color-map';
import { border, fillColor, pointRadius } from 'lib/deck/props/style';
import { selector } from 'recoil';
import {
  droughtOptionsVariableState,
  droughtRcpParamState,
  droughtRiskVariableState,
  droughtShowOptionsState,
  droughtShowRiskState,
} from 'state/drought/drought-parameters';
import { sectionVisibilityState } from 'state/sections';

export const droughtRegionsFieldSpecState = selector<FieldSpec>({
  key: 'droughtRegionsFieldSpecState',
  get: ({ get }) => {
    const field = get(droughtRiskVariableState);

    const rcp = get(droughtRcpParamState);

    return {
      fieldGroup: 'properties',
      field,
      fieldDimensions: DROUGHT_RISK_VARIABLES_WITH_RCP.includes(field)
        ? {
            rcp,
          }
        : {},
    };
  },
});

const droughtRiskColorSpecLookup: Record<DroughtRiskVariableType, ColorSpec> = {
  mean_monthly_water_stress_: VECTOR_COLOR_MAPS.droughtRiskWaterStress,
  epd: VECTOR_COLOR_MAPS.droughtRiskEpd,
  eael: VECTOR_COLOR_MAPS.droughtRiskEael,
};

export const droughtRegionsColorSpecState = selector<ColorSpec>({
  key: 'droughtRegionsColorSpecState',
  get: ({ get }) => {
    const field = get(droughtRiskVariableState);
    return droughtRiskColorSpecLookup[field];
  },
});

export const droughtRegionsLayerState = selector<ViewLayer>({
  key: 'droughtRegionsLayerState',
  get: ({ get }) => {
    const showDroughts = get(sectionVisibilityState('drought')) && get(droughtShowRiskState);

    if (!showDroughts) {
      return null;
    }

    const fieldSpec = get(droughtRegionsFieldSpecState);
    const colorSpec = get(droughtRegionsColorSpecState);

    const dataFn = getDroughtDataAccessor(fieldSpec);
    const colorFn = colorMap(colorSpec);

    return {
      id: 'drought_risk',
      group: null,
      interactionGroup: 'drought',
      dataAccessFn: getDroughtDataAccessor,
      dataFormatsFn: getDroughtRiskDataFormats,
      styleParams: {
        colorMap: {
          fieldSpec,
          colorSpec,
        },
      },

      fn: ({ deckProps, selection }) =>
        selectableMvtLayer(
          {
            selectionOptions: {
              selectedFeatureId: selection?.target?.feature.id,
            },
          },
          deckProps,
          {
            data: ASSETS_SOURCE.getDataUrl({ assetId: 'drought_combined' }),

            filled: true,
          },
          border([255, 255, 255]),
          fillColor(dataColorMap(dataFn, colorFn)),
        ),
    };
  },
});

export const droughtOptionsFieldSpecState = selector<FieldSpec>({
  key: 'droughtOptionsFieldSpecState',
  get: ({ get }) => {
    const field = get(droughtOptionsVariableState);

    const rcp = get(droughtRcpParamState);

    return {
      fieldGroup: 'properties',
      field,
      fieldDimensions: DROUGHT_OPTIONS_VARIABLES_WITH_RCP.includes(field)
        ? {
            rcp,
          }
        : {},
    };
  },
});

const droughtOptionsColorSpecLookup: Record<DroughtOptionsVariableType, ColorSpec> = {
  cost_jmd: VECTOR_COLOR_MAPS.droughtOptionsCost,
  population_protected: VECTOR_COLOR_MAPS.droughtOptionsPopulationProtected,
  net_present_value_benefit: VECTOR_COLOR_MAPS.droughtOptionsNPVBenefit,
  benefit_cost_ratio: VECTOR_COLOR_MAPS.droughtOptionsBenefitCost,
};

export const droughtOptionsColorSpecState = selector<ColorSpec>({
  key: 'droughtOptionsColorSpecState',
  get: ({ get }) => {
    const field = get(droughtOptionsVariableState);
    return droughtOptionsColorSpecLookup[field];
  },
});

export const droughtOptionsLayerState = selector<ViewLayer>({
  key: 'droughtOptionsLayerState',
  get: ({ get }) => {
    const showDroughts = get(sectionVisibilityState('drought')) && get(droughtShowOptionsState);

    if (!showDroughts) {
      return null;
    }

    const fieldSpec = get(droughtOptionsFieldSpecState);
    const colorSpec = get(droughtOptionsColorSpecState);

    const dataFn = getDroughtDataAccessor(fieldSpec);
    const colorFn = colorMap(colorSpec);

    return {
      id: 'drought_options',
      group: null,
      interactionGroup: 'drought',
      dataAccessFn: getDroughtDataAccessor,
      dataFormatsFn: getDroughtOptionsDataFormats,
      styleParams: {
        colorMap: {
          fieldSpec,
          colorSpec,
        },
      },
      fn: ({ deckProps, selection, zoom }) =>
        selectableMvtLayer(
          {
            selectionOptions: {
              selectedFeatureId: selection?.target?.feature.id,
            },
          },
          deckProps,
          {
            data: ASSETS_SOURCE.getDataUrl({ assetId: 'drought_options' }),

            filled: true,
          },
          pointRadius(zoom),
          border([255, 255, 255]),
          fillColor(dataColorMap(dataFn, colorFn)),
        ),
    };
  },
});
