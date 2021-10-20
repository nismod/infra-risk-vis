import { Box, Typography } from '@material-ui/core';
import { FC } from 'react';
import { LAYERS } from '../../config/layers';
// import { titleCase } from '../../helpers';
import { RasterHover, VectorHover } from '../DataMap';

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
      <span style={{color: logicalLayerSpec?.color ?? '#333'}}>■</span>&nbsp;
      <strong>
        {title} (ID: {f.properties.asset_id})
      </strong>
    </div>
  );
};

const RasterHoverDescription: FC<{ hoveredObject: RasterHover }> = ({ hoveredObject }) => {
  const { color } = hoveredObject;

  const sourceLogicalLayer = hoveredObject.logicalLayer;
  const logicalLayerSpec = LAYERS[sourceLogicalLayer];
  const title = logicalLayerSpec.label;

  return (
    <div>
      <span style={{color: `rgb(${color[0]},${color[1]},${color[2]})`}}>■</span>&nbsp;
      <strong>{title}</strong>
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
            <VectorHoverDescription hoveredObject={hv} key={hv.feature.id}/>
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
