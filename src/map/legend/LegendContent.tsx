import { Box, Typography } from '@material-ui/core';
import { useMemo } from 'react';
import * as d3Scale from 'd3-scale';
import * as d3Array from 'd3-array';

import { RASTER_COLOR_MAPS, VECTOR_COLOR_MAPS } from '../../config/color-maps';
import { DECK_LAYERS } from '../../config/deck-layers';
import { LAYERS } from '../../config/layers';
import { useRasterColorMapValues } from '../legend/use-color-map-values';

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

const RasterLegend = ({ deckLayerName, deckLayerParams }) => {
  const {
    sourceLogicalLayers: [logicalLayer],
    params,
  } = deckLayerParams;
  const logicalLayerConfig = LAYERS[logicalLayer];

  const { scheme, range } = RASTER_COLOR_MAPS[params.hazardType];

  const { error, loading, colorMapValues } = useRasterColorMapValues(scheme, range);

  return (
    <GradientLegend
      label={logicalLayerConfig.label}
      range={range}
      colorMapValues={!(error || loading) ? colorMapValues : null}
    />
  );
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

export const LegendContent = ({ deckLayersSpec, styleParams }) => {
  return (
    <>
      {/* <Typography variant="h6">Legend</Typography> */}
      {Object.entries(deckLayersSpec).map(([deckLayerName, params]) => {
        const deckLayerConfig = DECK_LAYERS[deckLayerName];

        if (deckLayerConfig.spatialType === 'raster') {
          return (
            <RasterLegend key={deckLayerName} deckLayerName={deckLayerName} deckLayerParams={params}></RasterLegend>
          );
        }

        return null;
      })}
      {styleParams.colorMap && <DamagesLegend styleParams={styleParams} />}
    </>
  );
};
