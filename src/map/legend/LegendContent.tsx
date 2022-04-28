import { FC, useCallback, useMemo } from 'react';

import { RASTER_COLOR_MAPS, VECTOR_COLOR_MAPS } from '../../config/color-maps';

import { useRasterColorMapValues } from '../legend/use-color-map-values';
import { HAZARDS_METADATA } from 'config/hazards/metadata';
import { ColorMap, FieldSpec, ViewLayer } from 'lib/data-map/view-layers';
import { useRecoilValue } from 'recoil';
import { viewLayersFlatState } from 'state/layers/view-layers-flat';
import { viewLayersParamsState } from 'state/layers/view-layers-params';
import { showPopulationState } from 'state/regions';
import { sectionVisibilityState } from 'state/sections';
import { colorScaleValues } from 'lib/color-map';
import { GradientLegend } from './GradientLegend';

const RasterLegend: FC<{ viewLayer: ViewLayer }> = ({ viewLayer }) => {
  const {
    params: { hazardType },
  } = viewLayer;
  const { label, dataUnit } = HAZARDS_METADATA[hazardType];
  const { scheme, range } = RASTER_COLOR_MAPS[hazardType];

  const { error, loading, colorMapValues } = useRasterColorMapValues(scheme, range);

  const getValueLabel = useCallback((value: number) => `${value.toLocaleString()} ${dataUnit}`, [dataUnit]);

  return (
    <GradientLegend
      label={label}
      range={range}
      colorMapValues={!(error || loading) ? colorMapValues : null}
      getValueLabel={getValueLabel}
    />
  );
};

interface LegendFormatParams {
  getLabel: (fieldSpec: FieldSpec) => string;
  getValueLabelFn: (fieldSpec: FieldSpec) => (value: number) => string;
}

const DAMAGES_LEGEND_PARAMS: LegendFormatParams = {
  getLabel: ({ field }) =>
    field.startsWith('ead') || field.startsWith('damages') ? 'Direct Damages' : 'Economic Losses',
  getValueLabelFn:
    ({ field }) =>
    (value: number) =>
      `$${value.toLocaleString()}`,
};

const POPULATION_LEGEND_PARAMS: LegendFormatParams = {
  getLabel: () => 'Population density',
  getValueLabelFn:
    ({ field }) =>
    (value: number) =>
      `${value.toLocaleString()}/km²`,
};

const ADAPTATION_LEGEND_PARAMS: LegendFormatParams = {
  getLabel: ({ field }) =>
    field === 'adaptation_cost'
      ? 'Adaptation Cost'
      : field === 'avoided_ead_mean'
      ? 'Avoided Damages'
      : field === 'avoided_eael_mean'
      ? 'Avoided Economic Losses'
      : 'Unknown',
  getValueLabelFn:
    ({ field }) =>
    (value: number) =>
      `$${value.toLocaleString()}`,
};

function getLegendFormatParams(viewLayer: ViewLayer, colorMap: ColorMap): LegendFormatParams {
  if (colorMap.fieldSpec.fieldGroup === 'damages_expected') {
    return DAMAGES_LEGEND_PARAMS;
  } else if (colorMap.fieldSpec.fieldGroup === 'adaptation') {
    return ADAPTATION_LEGEND_PARAMS;
  }
}

const VectorLegend: FC<{ colorMap: ColorMap; legendFormatParams: LegendFormatParams }> = ({
  colorMap,
  legendFormatParams,
}) => {
  const { colorSpec, fieldSpec } = colorMap;
  const colorMapValues = useMemo(() => colorScaleValues(colorSpec, 255), [colorSpec]);

  const { getLabel, getValueLabelFn } = legendFormatParams;

  const label = getLabel(fieldSpec);
  const getValueLabel = getValueLabelFn(fieldSpec);

  return (
    <GradientLegend
      label={label}
      range={colorSpec.range}
      colorMapValues={colorMapValues}
      getValueLabel={getValueLabel}
    />
  );
};

const populationColorMap: ColorMap = {
  colorSpec: VECTOR_COLOR_MAPS.population,
  fieldSpec: {} as FieldSpec,
};

export const LegendContent: FC<{}> = () => {
  const viewLayers = useRecoilValue(viewLayersFlatState);
  const viewLayersParams = useRecoilValue(viewLayersParamsState);
  const showRegions = useRecoilValue(sectionVisibilityState('regions'));
  const showPopulation = useRecoilValue(showPopulationState);

  const hazardViewLayers = [];

  let dataColorMaps: Record<string, [ColorMap, ViewLayer]> = {};

  viewLayers.forEach((viewLayer) => {
    if (viewLayer.spatialType === 'raster') {
      hazardViewLayers.push(viewLayer);
    } else {
      const { styleParams } = viewLayersParams[viewLayer.id];

      if (styleParams?.colorMap) {
        const { colorMap } = styleParams;
        const colorMapKey = `${colorMap.fieldSpec.fieldGroup}`;
        dataColorMaps[colorMapKey] ??= [colorMap, viewLayer];
      }
      // // save the first styleParams for damages
      // if (styleParams?.colorMap?.fieldSpec.fieldGroup === 'damages_expected' && !damagesColorMap) {
      //   damagesColorMap = styleParams.colorMap;
      // }
    }
  });

  return (
    <>
      {hazardViewLayers.map((viewLayer) => (
        <RasterLegend key={viewLayer.id} viewLayer={viewLayer} />
      ))}
      {Object.values(dataColorMaps).map(([colorMap, viewLayer]) => (
        <VectorLegend colorMap={colorMap} legendFormatParams={getLegendFormatParams(viewLayer, colorMap)} />
      ))}
      {/* {damagesColorMap && <VectorLegend colorMap={damagesColorMap} legendFormatParams={DAMAGES_LEGEND_PARAMS} />} */}
      {showRegions && showPopulation && (
        <VectorLegend colorMap={populationColorMap} legendFormatParams={POPULATION_LEGEND_PARAMS} />
      )}
    </>
  );
};
