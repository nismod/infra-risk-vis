import { FC } from 'react';
import { RegionsControl } from './RegionsControl';
import { SidebarPanel } from 'sidebar/SidebarPanel';
import { Box } from '@mui/material';
import { StyleSelection } from 'sidebar/StyleSelection';
import { regionsStyleState } from 'state/regions';
import { REGION_DEFAULT_STYLE, REGION_STYLES } from 'config/regions/styles';

export const RegionsSection: FC<{}> = () => {
  return (
    <SidebarPanel id="regions" title="Regions">
      <Box p={2}>
        <RegionsControl />
      </Box>
      <Box p={2} bgcolor="#eee">
        <StyleSelection state={regionsStyleState} options={REGION_STYLES} defaultValue={REGION_DEFAULT_STYLE} />
      </Box>
    </SidebarPanel>
  );
};
