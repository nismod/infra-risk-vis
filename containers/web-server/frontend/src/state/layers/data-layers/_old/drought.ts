import { selector } from 'recoil';

import { colorMap } from '@/lib/color-map';
import { ColorSpec, FieldSpec, ViewLayer } from '@/lib/data-map/view-layers';
import { selectableMvtLayer } from '@/lib/deck/layers/selectable-mvt-layer';
import { dataColorMap } from '@/lib/deck/props/color-map';
import { border, fillColor, pointRadius } from '@/lib/deck/props/style';

import { DROUGHT_OPTIONS_COLORMAPS, DROUGHT_RISK_COLORMAPS } from '@/config/_old/drought/colors';
import { getDroughtDataAccessor } from '@/config/_old/drought/data-access';
import { getDroughtOptionsDataFormats, getDroughtRiskDataFormats } from '@/config/_old/drought/data-formats';
import { DROUGHT_OPTIONS_VARIABLES_WITH_RCP, DROUGHT_RISK_VARIABLES_WITH_RCP } from '@/config/_old/drought/metadata';
import { SOURCES } from '@/config/sources';
import {
  droughtOptionsVariableState,
  droughtRcpParamState,
  droughtRiskVariableState,
  droughtShowOptionsState,
  droughtShowRiskState,
} from '@/state/data-selection/_old/drought/drought-parameters';
import { sectionVisibilityState } from '@/state/sections';

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

export const droughtRegionsColorSpecState = selector<ColorSpec>({
  key: 'droughtRegionsColorSpecState',
  get: ({ get }) => {
    const field = get(droughtRiskVariableState);
    return DROUGHT_RISK_COLORMAPS[field];
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
            data: SOURCES.vector.getUrl('drought_combined'),
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

export const droughtOptionsColorSpecState = selector<ColorSpec>({
  key: 'droughtOptionsColorSpecState',
  get: ({ get }) => {
    const field = get(droughtOptionsVariableState);
    return DROUGHT_OPTIONS_COLORMAPS[field];
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
            data: SOURCES.vector.getUrl('drought_options'),

            filled: true,
          },
          pointRadius(zoom),
          border([255, 255, 255]),
          fillColor(dataColorMap(dataFn, colorFn)),
        ),
    };
  },
});
