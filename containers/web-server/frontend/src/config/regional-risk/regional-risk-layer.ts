import * as d3Scale from 'd3-scale';
import * as d3ScaleChromatic from 'd3-scale-chromatic';
import React from 'react';

import { colorMap } from '@/lib/color-map';
import { InteractionTarget, VectorTarget } from '@/lib/data-map/interactions/use-interactions';
import { ColorSpec, FieldSpec, ViewLayer } from '@/lib/data-map/view-layers';
import { selectableMvtLayer } from '@/lib/deck/layers/selectable-mvt-layer';
import { dataColorMap } from '@/lib/deck/props/color-map';
import { featureProperty } from '@/lib/deck/props/data-source';
import { border, fillColor } from '@/lib/deck/props/style';

import { SimpleAssetDetails } from '@/details/features/asset-details';
import { VectorHoverDescription } from '@/map/tooltip/VectorHoverDescription';

import { SOURCES } from '../sources';
import { getRegionalExposureDataFormats } from './data-formats';
import { RegionalExposureDetails } from './details';
import { RegionalExposureVariableType } from './metadata';

const rexpColorLookup: Record<RegionalExposureVariableType, ColorSpec> = {
  'pop_exposed_seismic_threshold0.1g': {
    scale: d3Scale.scaleSequential,
    scheme: d3ScaleChromatic.interpolateBlues,
    range: [0, 100_000_000],
    empty: '#ccc',
  },
  'pop_exposed_seismic_threshold0.2g': {
    scale: d3Scale.scaleSequential,
    scheme: d3ScaleChromatic.interpolateBlues,
    range: [0, 100_000_000],
    empty: '#ccc',
  },
  pop_exposed_river_historical_WATCH_1980_thresholdNone: {
    scale: d3Scale.scaleSequential,
    scheme: d3ScaleChromatic.interpolateBlues,
    range: [0, 100_000_000],
    empty: '#ccc',
  },
  'pop_exposed_river_rcp4p5_MIROC-ESM-CHEM_2050_thresholdNone': {
    scale: d3Scale.scaleSequential,
    scheme: d3ScaleChromatic.interpolateBlues,
    range: [0, 100_000_000],
    empty: '#ccc',
  },
};

export function regionalExposureLayer(variable: RegionalExposureVariableType): ViewLayer {
  const fieldSpec: FieldSpec = {
    fieldGroup: 'properties',
    field: variable,
  };

  const colorSpec = rexpColorLookup[variable];

  const id = `adm0_exposure`;

  return {
    id,
    interactionGroup: 'rexp',
    spatialType: 'vector',
    params: {
      variable,
    },
    styleParams: {
      colorMap: {
        fieldSpec,
        colorSpec,
      },
    },
    fn: ({ deckProps, zoom, selection }) => {
      const dataStyleColor = dataColorMap(featureProperty(variable), colorMap(colorSpec));

      return selectableMvtLayer(
        {
          selectionOptions: {
            selectedFeatureId: selection?.target.feature.id,
            selectionFillColor: [0, 0, 0, 0],
            selectionLineColor: [0, 255, 255, 255],
          },
        },
        deckProps,
        {
          data: SOURCES.vector.getUrl(id),
        },
        border([100, 100, 100]),
        fillColor(dataStyleColor),
        {
          highlightColor: [255, 255, 255, 100],
        },
      );
    },
    dataFormatsFn: getRegionalExposureDataFormats,
    dataAccessFn: ({ field }) => featureProperty(field),
    renderDetails(selection: InteractionTarget<VectorTarget>) {
      const feature = selection.target.feature;

      return React.createElement(SimpleAssetDetails, {
        feature: feature,
        label: 'Regional Exposure',
        detailsComponent: RegionalExposureDetails,
      });
    },
    renderTooltip: (hover: InteractionTarget<VectorTarget>) => {
      return React.createElement(VectorHoverDescription, {
        hoveredObject: hover,
        label: 'Regional Exposure',
        color: '#83B4FF',
        idValue: hover.target.feature.properties.ISO_A3,
      });
    },
  };
}
