import { useRecoilState } from 'recoil';

import { ParamChecklist } from '@/lib/controls/params/ParamChecklist';

import { INDUSTRY_METADATA } from '@/config/industry/industry-view-layer';
import { LayerLabel, LayerLabelShapeType } from '@/sidebar/ui/LayerLabel';
import { IndustryType, industrySelectionState } from '@/state/data-selection/industry';

export const IndustryControl = () => {
  const [checkboxState, setCheckboxState] = useRecoilState(industrySelectionState);

  return (
    <ParamChecklist<IndustryType>
      title="Industry types"
      options={Object.keys(checkboxState) as IndustryType[]}
      checklistState={checkboxState}
      onChecklistState={setCheckboxState}
      renderLabel={(key) => {
        const { color, type, label } = INDUSTRY_METADATA[key];
        return <LayerLabel color={color} type={type as LayerLabelShapeType} label={label} />;
      }}
    />
  );
};
