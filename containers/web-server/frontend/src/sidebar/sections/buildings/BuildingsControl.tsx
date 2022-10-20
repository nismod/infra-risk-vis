import { useRecoilState } from 'recoil';

import { ParamChecklist } from '@/lib/controls/params/ParamChecklist';

import { NETWORKS_METADATA } from '@/config/networks/metadata';
import { LayerLabel } from '@/sidebar/ui/LayerLabel';
import { buildingSelectionState } from '@/state/buildings';

export const BuildingsControl = () => {
  const [checkboxState, setCheckboxState] = useRecoilState(buildingSelectionState);

  return (
    <ParamChecklist
      title="Building types"
      options={Object.keys(checkboxState)}
      checklistState={checkboxState}
      onChecklistState={setCheckboxState}
      renderLabel={(key) => <LayerLabel {...NETWORKS_METADATA[key]} />}
    />
  );
};
