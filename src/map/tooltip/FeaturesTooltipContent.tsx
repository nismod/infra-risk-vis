import { Feature } from 'geojson';
import { FC } from 'react';

import { titleCase } from '../../helpers';

export const FeaturesTooltipContent: FC<{
  features: Feature[];
}> = ({ features }) => {
  const entries: object = {};
  for (const f of features) {
    const sourceLayer = f.properties.layerName;
    let title = titleCase(
      sourceLayer.replace(/_/g, ' ').replace('edges', '').replace('nodes', '').replace('elec', 'electricity'),
    );
    let subtitle = f.properties.road_type ? '(' + f.properties.road_type + ')' : '';

    if (!entries[sourceLayer]) {
      entries[sourceLayer] = { title, subtitle };
    }
  }

  return features.length ? (
    <>
      {Object.values(entries).map((entry, i) => {
        return (
          <div key={i}>
            <strong>
              {entry.title} {entry.subtitle}
            </strong>
            {entry.detail}
          </div>
        );
      })}
    </>
  ) : null;
};
