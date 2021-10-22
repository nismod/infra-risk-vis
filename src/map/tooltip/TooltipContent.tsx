import { Box, Typography } from '@material-ui/core';
import { FC, useMemo } from 'react';
import { RASTER_COLOR_MAPS } from '../../config/color-maps';
import { DECK_LAYERS } from '../../config/deck-layers';
import { LAYERS } from '../../config/layers';
// import { titleCase } from '../../helpers';
import { RasterHover, VectorHover } from '../DataMap';
import { useRasterColorMapValues } from '../legend/use-color-map-values';

const VectorHoverDescription: FC<{ hoveredObject: VectorHover }> = ({ hoveredObject }) => {
  const f = hoveredObject.feature;
  // const sourceDeckLayer = hoveredObject.deckLayer;
  const sourceLogicalLayer = hoveredObject.logicalLayer;
  const logicalLayerSpec = LAYERS[sourceLogicalLayer];
  const title = logicalLayerSpec.label;

  // let title = titleCase(
  //   sourceLayer.replace(/_/g, ' ').replace('edges', '').replace('nodes', '').replace('elec', 'electricity'),
  // );
  // let subtitle = f.properties.road_type ? '(' + f.properties.road_type + ')' : '';

  // if (!entries[sourceLayer]) {
  //   entries[sourceLayer] = { title, subtitle };
  // }

  return (
    <div>
      <span style={{ color: logicalLayerSpec?.color ?? '#333' }}>■</span>&nbsp;
      <strong>
        {title} (ID: {f.properties.asset_id})
      </strong>
    </div>
  );
};

function useRasterColorMapLookup(colorMapValues) {
  return useMemo(
    () => colorMapValues && Object.fromEntries(colorMapValues.map(({ value, color }) => [color, value])),
    [colorMapValues],
  );
}

const RasterHoverDescription: FC<{ hoveredObject: RasterHover }> = ({ hoveredObject }) => {
  const { color } = hoveredObject;

  const { logicalLayer, deckLayer } = hoveredObject;
  const { label, dataUnit } = LAYERS[logicalLayer];
  const {
    dataParams: { hazardType },
  } = DECK_LAYERS[deckLayer];
  const { scheme, range } = RASTER_COLOR_MAPS[hazardType];

  const title = `${label} (${dataUnit})`;

  const { colorMapValues } = useRasterColorMapValues(scheme, range);
  const rasterValueLookup = useRasterColorMapLookup(colorMapValues);

  const colorString = `rgb(${color[0]},${color[1]},${color[2]})`;
  return (
    <div>
      <span style={{ color: colorString }}>■</span>&nbsp;
      <strong>{title}</strong>
      {rasterValueLookup?.[colorString] && <span>: {rasterValueLookup[colorString].toFixed(2)}</span>}
    </div>
  );
};

export const TooltipContent: FC<{ hoveredVectors: VectorHover[]; hoveredRasters: RasterHover[] }> = ({
  hoveredVectors,
  hoveredRasters,
}) => {
  return (
    <>
      {hoveredVectors.length ? (
        <Box mb={2}>
          <Typography>Asset</Typography>
          {hoveredVectors.map((hv) => (
            <VectorHoverDescription hoveredObject={hv} key={hv.feature.id} />
          ))}
        </Box>
      ) : null}
      {hoveredRasters.length ? (
        <Box>
          <Typography>Hazards</Typography>
          {hoveredRasters.map((hr) => (
            <RasterHoverDescription hoveredObject={hr} key={`${hr.deckLayer}-${hr.color}`} />
          ))}
        </Box>
      ) : null}
    </>
  );
};
