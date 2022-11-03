import { useRecoilState } from 'recoil';

import { ParamChecklist } from '@/lib/controls/params/ParamChecklist';
import { ValueLabel } from '@/lib/controls/params/value-label';

import { NATURE_RASTER_VALUE_LABELS } from '@/config/natural-assets/metadata';
import { NaturalAssetType, naturalAssetsSelectionState } from '@/state/data-selection/natural-assets';

const NATURAL_ASSET_VALUE_LABELS = NATURE_RASTER_VALUE_LABELS.filter(
  (x) => x.value === 'organic_carbon',
) as ValueLabel<NaturalAssetType>[];

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
