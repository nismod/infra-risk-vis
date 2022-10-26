import { useRecoilState } from 'recoil';

import { ParamChecklist } from '@/lib/controls/params/ParamChecklist';

import { NETWORKS_METADATA } from '@/config/networks/metadata';
import { LayerLabel, LayerLabelShapeType } from '@/sidebar/ui/LayerLabel';
import { IndustryType, industrySelectionState } from '@/state/industry';

export const IndustryControl = () => {
  const [checkboxState, setCheckboxState] = useRecoilState(industrySelectionState);

  return (
    <ParamChecklist<IndustryType>
      title="Industry types"
      options={Object.keys(checkboxState) as IndustryType[]}
      checklistState={checkboxState}
      onChecklistState={setCheckboxState}
      renderLabel={(key) => {
        const { color, type, label } = NETWORKS_METADATA[key];
        return <LayerLabel color={color} type={type as LayerLabelShapeType} label={label} />;
      }}
    />
  );
};
