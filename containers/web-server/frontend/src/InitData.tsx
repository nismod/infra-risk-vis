import { Suspense } from 'react';

import { HAZARDS_UI_ORDER } from './config/hazards/metadata';
import { InitHazardData } from './sidebar/sections/hazards/HazardsControl';

export const InitData = () => {
  return (
    <Suspense fallback={null}>
      {HAZARDS_UI_ORDER.map((hazard) => (
        <InitHazardData key={hazard} type={hazard} />
      ))}
    </Suspense>
  );
};
