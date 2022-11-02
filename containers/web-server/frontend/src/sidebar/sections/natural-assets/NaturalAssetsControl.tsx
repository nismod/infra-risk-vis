import { useRecoilState } from 'recoil';

import { ParamChecklist } from '@/lib/controls/params/ParamChecklist';

import { NATURAL_ASSET_VALUE_LABELS, NaturalAssetType } from '@/config/natural-assets/metadata';
import { naturalAssetsSelectionState } from '@/state/data-selection/natural-assets';

export const NaturalAssetsControl = () => {
  const [checklistState, setChecklistState] = useRecoilState(naturalAssetsSelectionState);
  return (
    <ParamChecklist<NaturalAssetType>
      title={null}
      options={NATURAL_ASSET_VALUE_LABELS}
      checklistState={checklistState}
      onChecklistState={setChecklistState}
      renderLabel={(key, label) => <>{label}</>}
      showAllNone={false}
    />
  );
};
