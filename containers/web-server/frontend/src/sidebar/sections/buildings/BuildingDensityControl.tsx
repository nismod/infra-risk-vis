import { useRecoilState } from 'recoil';

import { ParamDropdown } from '@/lib/controls/ParamDropdown';

import {
  BUILDING_DENSITY_TYPE_LABELS,
  BuildingDensityType,
} from '@/config/building-density/metadata';
import { buildingDensityTypeState } from '@/state/data-selection/building-density';

export const BuildingDensityControl = () => {
  const [type, setType] = useRecoilState(buildingDensityTypeState);
  return (
    <ParamDropdown<BuildingDensityType>
      value={type}
      onChange={setType}
      options={BUILDING_DENSITY_TYPE_LABELS}
    />
  );
};
