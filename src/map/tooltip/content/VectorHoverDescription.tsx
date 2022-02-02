import { FC } from 'react';

import { LAYERS } from '../../../config/layers';
import { VectorHover } from '../../DataMap';

export const VectorHoverDescription: FC<{ hoveredObject: VectorHover }> = ({ hoveredObject }) => {
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