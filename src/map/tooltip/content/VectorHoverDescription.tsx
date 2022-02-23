import { NETWORKS_METADATA } from 'config/networks/metadata';
import { InteractionTarget, VectorTarget } from 'lib/data-map/interactions/use-interactions';
import { FC } from 'react';

export const VectorHoverDescription: FC<{ hoveredObject: InteractionTarget<VectorTarget> }> = ({ hoveredObject }) => {
  const f = hoveredObject.target.feature;
  // const sourceDeckLayer = hoveredObject.deckLayer;
  const { label: title, color } = NETWORKS_METADATA[hoveredObject.viewLayer.id];

  // let title = titleCase(
  //   sourceLayer.replace(/_/g, ' ').replace('edges', '').replace('nodes', '').replace('elec', 'electricity'),
  // );
  // let subtitle = f.properties.road_type ? '(' + f.properties.road_type + ')' : '';

  // if (!entries[sourceLayer]) {
  //   entries[sourceLayer] = { title, subtitle };
  // }

  return (
    <div>
      <span style={{ color: color ?? '#333' }}>■</span>&nbsp;
      <strong>
        {title} (ID: {f.properties.asset_id})
      </strong>
    </div>
  );
};
