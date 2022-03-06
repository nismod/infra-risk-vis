import { Box } from '@mui/material';
import { StyleSelection } from 'sidebar/StyleSelection';
import { buildingsStyleState } from 'state/buildings';
import { SidebarPanel } from '../SidebarPanel';

export const BuildingsSection = () => {
  return (
    <SidebarPanel id="buildings" title="Buildings">
      <Box p={2} bgcolor="#eee">
        <StyleSelection
          state={buildingsStyleState}
          options={[
            {
              label: 'Building type',
              value: 'type',
            },
          ]}
          defaultValue="type"
        />
      </Box>
    </SidebarPanel>
  );
};
