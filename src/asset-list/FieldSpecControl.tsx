import { Box, Typography } from '@mui/material';
import { ParamDropdown } from 'lib/controls/ParamDropdown';
import { useCallback } from 'react';

export const FieldSpecControl = ({ fieldSpec, onFieldSpec }) => {
  const { field, fieldParams } = fieldSpec;

  const handleField = useCallback((value) => onFieldSpec({ field: value, fieldParams }), [onFieldSpec, fieldParams]);
  const handleFieldParams = useCallback((value) => onFieldSpec({ field, fieldParams: value }), [onFieldSpec, field]);

  return (
    <>
      <Box mt={1}>
        <ParamDropdown
          title="Variable"
          value={field}
          onChange={handleField}
          options={[{ value: 'damages', label: 'Damages' }]}
        />
      </Box>
      <Box mt={1}>
        <Typography variant="subtitle1">Parameters</Typography>
        <FieldParamsControl field={field} fieldParams={fieldParams} onFieldParams={handleFieldParams} />
      </Box>
    </>
  );
};

function FieldParamsControl({ field, fieldParams, onFieldParams }) {
  const handleParamChange = useCallback(
    (param, value) => {
      onFieldParams({ ...fieldParams, [param]: value });
    },
    [onFieldParams, fieldParams],
  );

  return (
    <>
      <Box>
        <ParamDropdown
          title="Hazard"
          value={fieldParams.hazard}
          onChange={(value) => handleParamChange('hazard', value)}
          options={[
            { value: 'total-damages', label: 'All Hazards' },
            { value: 'fluvial', label: 'River Flooding' },
            { value: 'surface', label: 'Surface Flooding' },
            { value: 'coastal', label: 'Coastal Flooding' },
            { value: 'cyclone', label: 'Cyclones' },
          ]}
        />
        <ParamDropdown
          title="RCP"
          value={fieldParams.rcp}
          onChange={(value) => handleParamChange('rcp', value)}
          options={[{ value: 'baseline', label: 'Baseline' }, '2.6', '4.5', '8.5']}
        />
        <ParamDropdown
          title="Epoch"
          value={fieldParams.epoch}
          onChange={(value) => handleParamChange('epoch', value)}
          options={[2010, 2050, 2100]}
        />
      </Box>
    </>
  );
}
