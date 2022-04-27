import { Box, Typography } from '@mui/material';
import { FC, useCallback, useMemo } from 'react';
import * as d3Scale from 'd3-scale';
import * as d3Array from 'd3-array';

import { RASTER_COLOR_MAPS, VECTOR_COLOR_MAPS } from '../../config/color-maps';

import { useRasterColorMapValues } from '../legend/use-color-map-values';
import { HAZARDS_METADATA } from 'config/hazards/metadata';
import { ViewLayer } from 'lib/data-map/view-layers';
import { useRecoilValue } from 'recoil';
import { viewLayersFlatState } from 'state/layers/view-layers-flat';
import { viewLayersParamsState } from 'state/layers/view-layers-params';
import { showPopulationState } from 'state/regions';
import { sectionVisibilityState } from 'state/sections';

const legendHeight = 10;

const LegendGradient: FC<{
  colorMapValues: any[];
  getValueLabel: (value: number) => string;
}> = ({ colorMapValues, getValueLabel }) => {
  return (
    <>
      {colorMapValues.map(({ color, value }, i) => (
        <Box key={i} height={legendHeight} width={1} bgcolor={color} title={getValueLabel(value)} />
      ))}
    </>
  );
};

const GradientLegend = ({ label, range, colorMapValues, getValueLabel }) => (
  <Box mb={2}>
    <Typography>{label}</Typography>
    <Box
      height={legendHeight + 2}
      width={255}
      bgcolor="#ccc"
      display="flex"
      flexDirection="row"
      border="1px solid gray"
    >
      {colorMapValues && <LegendGradient colorMapValues={colorMapValues} getValueLabel={getValueLabel} />}
    </Box>
    <Box height={10} position="relative">
      {colorMapValues && (
        <>
          <Box position="absolute" left={0}>
            <Typography>{getValueLabel(range[0])}</Typography>
          </Box>
          <Box position="absolute" right={0}>
            <Typography>{getValueLabel(range[1])}</Typography>
          </Box>
        </>
      )}
    </Box>
  </Box>
);

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

const DamagesLegend = ({ styleParams }) => {
  const {
    colorMap: { colorScheme, colorField },
  } = styleParams;

  const { scale, range } = VECTOR_COLOR_MAPS[colorScheme];
  const [rangeMin, rangeMax] = range;

  const { field } = colorField;

  const isDirect = field.startsWith('ead') || field.startsWith('damages');

  const colorMapValues = useMemo(() => {
    const scaleFn = d3Scale.scaleSequential([rangeMin, rangeMax], scale);

    return d3Array.ticks(rangeMin, rangeMax, 255).map((x) => ({ value: x, color: scaleFn(x) }));
  }, [scale, rangeMin, rangeMax]);

  const getValueLabel = useCallback((value: number) => `${value.toLocaleString()}$`, []);

  // const { error, loading, colorMapValues } = useVectorColorMapValues(scheme, range);
  const label = isDirect ? 'Direct Damages' : 'Economic Losses';

  return <GradientLegend label={label} range={range} colorMapValues={colorMapValues} getValueLabel={getValueLabel} />;
};

const PopulationLegend = () => {
  const { scale, range } = VECTOR_COLOR_MAPS['population'];
  const [rangeMin, rangeMax] = range;

  const colorMapValues = useMemo(() => {
    const scaleFn = d3Scale.scaleSequential([rangeMin, rangeMax], scale);

    return d3Array.ticks(rangeMin, rangeMax, 255).map((x) => ({ value: x, color: scaleFn(x) }));
  }, [scale, rangeMin, rangeMax]);

  const getValueLabel = useCallback((value: number) => `${value.toLocaleString()}/km²`, []);

  // const { error, loading, colorMapValues } = useVectorColorMapValues(scheme, range);

  return (
    <GradientLegend
      label="Population density"
      range={range}
      colorMapValues={colorMapValues}
      getValueLabel={getValueLabel}
    />
  );
};

export const LegendContent: FC<{}> = () => {
  const viewLayers = useRecoilValue(viewLayersFlatState);
  const viewLayersParams = useRecoilValue(viewLayersParamsState);
  const showRegions = useRecoilValue(sectionVisibilityState('regions'));
  const showPopulation = useRecoilValue(showPopulationState);

  const hazardViewLayers = [];
  let damageStyleParams = null;

  viewLayers.forEach((viewLayer) => {
    if (viewLayer.spatialType === 'raster') {
      hazardViewLayers.push(viewLayer);
    } else {
      const { styleParams } = viewLayersParams[viewLayer.id];

      // save the first styleParams for damages
      if (styleParams?.colorMap?.colorScheme === 'damages' && !damageStyleParams) {
        damageStyleParams = styleParams;
      }
    }
  });

  return (
    <>
      {hazardViewLayers.map((viewLayer) =>
        viewLayer.spatialType === 'raster' ? <RasterLegend key={viewLayer.id} viewLayer={viewLayer} /> : null,
      )}
      {damageStyleParams && <DamagesLegend styleParams={damageStyleParams} />}
      {showRegions && showPopulation && <PopulationLegend />}
    </>
  );
};
