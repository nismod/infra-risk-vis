import { FC } from 'react';
import { LAYERS } from '../../config/layers';
import { titleCase } from '../../helpers';
import { HoveredObject, RasterHover, VectorHover } from '../DataMap';

const VectorHoverDescription: FC<{ hoveredObject: VectorHover }> = ({ hoveredObject }) => {
  const f = hoveredObject.feature;
  const sourceLayer = hoveredObject.deckLayer;
  const title = f.properties.layerName;

  // let title = titleCase(
  //   sourceLayer.replace(/_/g, ' ').replace('edges', '').replace('nodes', '').replace('elec', 'electricity'),
  // );
  // let subtitle = f.properties.road_type ? '(' + f.properties.road_type + ')' : '';

  // if (!entries[sourceLayer]) {
  //   entries[sourceLayer] = { title, subtitle };
  // }

  const logicalLayerSpec = LAYERS[sourceLayer];

  return (
    <div key={`${title}-${f.id}`}>
      <strong style={{ color: logicalLayerSpec?.color ?? '#333' }}>
        {title} ({f.id})
      </strong>
    </div>
  );
};

const RasterHoverDescription: FC<{ hoveredObject: RasterHover }> = ({ hoveredObject }) => {
  const { color, deckLayer } = hoveredObject;

  return (
    <div key={`${deckLayer}-${color}`} style={{ backgroundColor: `rgb(${color[0]},${color[1]},${color[2]})` }}>
      <strong>{deckLayer}</strong>
    </div>
  );
};

export const TooltipContent: FC<{ hoveredObjects: HoveredObject[] }> = ({ hoveredObjects }) => {
  return (
    <>
      {hoveredObjects.map((ho) =>
        ho.type === 'vector' ? (
          <VectorHoverDescription hoveredObject={ho} />
        ) : (
          <RasterHoverDescription hoveredObject={ho} />
        ),
      )}
    </>
  );
};
