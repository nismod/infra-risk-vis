import { useRecoilState } from 'recoil';

import { ParamChecklist } from '@/lib/controls/params/ParamChecklist';

import { BUILDINGS_METADATA } from '@/config/_old/buildings/metadata';
import { LayerLabel, LayerLabelShapeType } from '@/sidebar/ui/LayerLabel';
import { BuildingType, buildingSelectionState } from '@/state/data-selection/_old/buildings';

export const BuildingsControl = () => {
  const [checkboxState, setCheckboxState] = useRecoilState(buildingSelectionState);

  return (
    <ParamChecklist<BuildingType>
      title="Building types"
      options={Object.keys(checkboxState) as BuildingType[]}
      checklistState={checkboxState}
      onChecklistState={setCheckboxState}
      renderLabel={(key) => {
        const { label, type, color } = BUILDINGS_METADATA[key];
        return <LayerLabel {...{ label, type, color }} />;
      }}
    />
  );
};
