import { Box, Typography } from '@material-ui/core';

import { RASTER_COLOR_MAPS } from '../../config/color-maps';
import { DECK_LAYERS } from '../../config/deck-layers';
import { LAYERS } from '../../config/layers';
import { useColorMapValues } from '../legend/use-color-map-values';

const legendHeight = 10;

const RasterLegendGradient = ({ colorMapValues }) => {
  return (
    <>
      {colorMapValues.map(({ rgba: [r, g, b], value }) => (
        <Box
          // key={`${value}-${r}-${g}-${b}`}
          height={legendHeight}
          width={1}
          bgcolor={`rgb(${r},${g},${b})`}
          title={value.toFixed(3)}
        ></Box>
      ))}
    </>
  );
};

const RasterLegend = ({ deckLayerName, deckLayerParams }) => {
  const {
    sourceLogicalLayers: [logicalLayer],
    params,
  } = deckLayerParams;
  const logicalLayerConfig = LAYERS[logicalLayer];

  const { scheme, range } = RASTER_COLOR_MAPS[params.hazardType];

  const { error, loading, colorMapValues } = useColorMapValues(scheme, range);

  return (
    <Box mb={2}>
      <Typography>{logicalLayerConfig.label}</Typography>
      <Box
        height={legendHeight + 2}
        width={255}
        bgcolor="#ccc"
        display="flex"
        flexDirection="row"
        border="1px solid gray"
      >
        {!loading && !error && <RasterLegendGradient colorMapValues={colorMapValues} />}
      </Box>
      <Box height={10} position="relative">
        {!loading && !error && (
          <>
            <Box position="absolute" left={0}>
              <Typography>{range[0]}</Typography>
            </Box>
            <Box position="absolute" right={0}>
              <Typography>{range[1]}</Typography>
            </Box>
          </>
        )}
      </Box>
    </Box>
  );
};

export const LegendContent = ({ deckLayersSpec }) => {
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
    </>
  );
};
