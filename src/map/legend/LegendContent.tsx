import { Box, Typography } from '@mui/material';
import { FC, useMemo } from 'react';
import * as d3Scale from 'd3-scale';
import * as d3Array from 'd3-array';

import { RASTER_COLOR_MAPS, VECTOR_COLOR_MAPS } from '../../config/color-maps';

import { useRasterColorMapValues } from '../legend/use-color-map-values';
import { HAZARDS_METADATA } from 'config/hazards/metadata';
import { ViewLayer, ViewLayerParams } from 'lib/data-map/view-layers';

const legendHeight = 10;

const LegendGradient = ({ colorMapValues }) => {
  return (
    <>
      {colorMapValues.map(({ color, value }) => (
        <Box
          // key={`${value}-${r}-${g}-${b}`}
          height={legendHeight}
          width={1}
          bgcolor={color}
          title={value.toFixed(3)}
        ></Box>
      ))}
    </>
  );
};

const GradientLegend = ({ label, range, colorMapValues }) => (
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
      {colorMapValues && <LegendGradient colorMapValues={colorMapValues} />}
    </Box>
    <Box height={10} position="relative">
      {colorMapValues && (
        <>
          <Box position="absolute" left={0}>
            <Typography>{range[0].toLocaleString()}</Typography>
          </Box>
          <Box position="absolute" right={0}>
            <Typography>{range[1].toLocaleString()}</Typography>
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
  const { label } = HAZARDS_METADATA[hazardType];
  const { scheme, range } = RASTER_COLOR_MAPS[hazardType];

  const { error, loading, colorMapValues } = useRasterColorMapValues(scheme, range);

  return <GradientLegend label={label} range={range} colorMapValues={!(error || loading) ? colorMapValues : null} />;
};

const DamagesLegend = ({ styleParams }) => {
  const {
    colorMap: { colorScheme },
  } = styleParams;

  const { scale, range } = VECTOR_COLOR_MAPS[colorScheme];
  const [rangeMin, rangeMax] = range;

  const colorMapValues = useMemo(() => {
    const scaleFn = d3Scale.scaleSequential([rangeMin, rangeMax], scale);

    return d3Array.ticks(rangeMin, rangeMax, 255).map((x) => ({ value: x, color: scaleFn(x) }));
  }, [scale, rangeMin, rangeMax]);

  // const { error, loading, colorMapValues } = useVectorColorMapValues(scheme, range);

  return <GradientLegend label="Direct Damages" range={range} colorMapValues={colorMapValues} />;
};

export const LegendContent: FC<{ viewLayers: ViewLayer[]; viewLayersParams: Record<string, ViewLayerParams> }> = ({
  viewLayers,
  viewLayersParams,
}) => {
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
    </>
  );
};
