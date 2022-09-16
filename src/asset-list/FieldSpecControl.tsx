import { Box, Typography } from '@mui/material';
import produce from 'immer';
import { ParamDropdown } from 'lib/controls/ParamDropdown';
import { FieldSpec } from 'lib/data-map/view-layers';
import { FC } from 'react';

export const FieldSpecControl: FC<{
  fieldSpec: FieldSpec;
  onFieldSpec: (fieldSpec: FieldSpec) => void;
}> = ({ fieldSpec, onFieldSpec }) => {
  const { fieldGroup, fieldDimensions } = fieldSpec;

  return (
    <>
      <Box mt={1}>
        <ParamDropdown
          title="Variable Group"
          value={fieldGroup}
          onChange={(value) =>
            onFieldSpec(
              produce(fieldSpec, (draft) => {
                draft.fieldGroup = value;
              }),
            )
          }
          options={[
            { value: 'damages_expected', label: 'Expected Damages' },
            { value: 'adaptation', label: 'Adaptation Options' },
          ]}
        />
      </Box>
      <Box mt={1}>
        <Typography variant="h6">Dimensions</Typography>
        <Box>
          <ParamDropdown
            title="Hazard"
            value={fieldDimensions.hazard}
            onChange={(value) =>
              onFieldSpec(
                produce(fieldSpec, (draft) => {
                  draft.fieldDimensions.hazard = value;
                }),
              )
            }
            options={[
              { value: 'all', label: 'All Hazards' },
              { value: 'river', label: 'River Flooding' },
              { value: 'coastal', label: 'Coastal Flooding' },
            ]}
          />
          <ParamDropdown
            title="RCP"
            value={fieldDimensions.rcp}
            onChange={(value) =>
              onFieldSpec(
                produce(fieldSpec, (draft) => {
                  draft.fieldDimensions.rcp = value;
                }),
              )
            }
            options={[{ value: 'baseline', label: 'Baseline' }, '4.5', '8.5']}
          />
          <ParamDropdown
            title="Epoch"
            value={fieldDimensions.epoch}
            onChange={(value) =>
              onFieldSpec(
                produce(fieldSpec, (draft) => {
                  draft.fieldDimensions.epoch = value;
                }),
              )
            }
            options={[2010, 2030, 2050, 2080]}
          />
        </Box>
      </Box>
    </>
  );
};
