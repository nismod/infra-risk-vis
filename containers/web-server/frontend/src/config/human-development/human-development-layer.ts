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
import { toLabelLookup } from '@/lib/helpers';

import { SimpleAssetDetails } from '@/details/features/asset-details';

import { SOURCES } from '../sources';
import { getHumanDevelopmentDataFormats } from './data-formats';
import { HDI_REGION_LEVEL_DETAILS } from './details';
import { HDI_REGION_LEVEL_LABELS, HdiRegionLevel, HdiVariableType } from './metadata';

const hdiColorLookup: Record<HdiVariableType, ColorSpec> = {
  subnational_hdi: {
    scale: d3Scale.scaleSequential,
    scheme: d3ScaleChromatic.interpolateBlues,
    range: [0, 1],
    empty: '#ccc',
  },
  educational_index: {
    scale: d3Scale.scaleSequential,
    scheme: d3ScaleChromatic.interpolatePurples,
    range: [0, 1],
    empty: '#ccc',
  },
  health_index: {
    scale: d3Scale.scaleSequential,
    scheme: d3ScaleChromatic.interpolateGreens,
    range: [0, 1],
    empty: '#ccc',
  },
  income_index: {
    scale: d3Scale.scaleSequential,
    scheme: d3ScaleChromatic.interpolateOranges,
    range: [0, 1],
    empty: '#ccc',
  },
};

const hdiRegionLabelLookup = toLabelLookup(HDI_REGION_LEVEL_LABELS);

export function humanDevelopmentLayer(regionLevel: HdiRegionLevel, variable: HdiVariableType): ViewLayer {
  const fieldSpec: FieldSpec = {
    fieldGroup: 'properties',
    field: variable,
  };

  const colorSpec = hdiColorLookup[variable];
  const regionLevelLabel = hdiRegionLabelLookup[regionLevel];

  const id = `hdi_${regionLevel}`;

  return {
    id,
    interactionGroup: 'hdi',
    spatialType: 'vector',
    params: {
      regionLevel,
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
        border([40, 40, 40, 255]),
        fillColor(dataStyleColor),
        {
          highlightColor: [255, 255, 255, 100],
        },
      );
    },
    dataFormatsFn: getHumanDevelopmentDataFormats,
    dataAccessFn: ({ field }) => featureProperty(field),
    renderDetails(selection: InteractionTarget<VectorTarget>) {
      const feature = selection.target.feature;
      const detailsComponent = HDI_REGION_LEVEL_DETAILS[regionLevel];

      return React.createElement(SimpleAssetDetails, {
        detailsComponent,
        feature: feature,
        label: `Human Development (${regionLevelLabel})`,
      });
    },
  };
}
