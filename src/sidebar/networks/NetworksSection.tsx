import { FC } from 'react';
import { NetworkControl } from './NetworkControl';
import { SidebarPanel } from 'sidebar/SidebarPanel';
import { Box } from '@mui/material';
import { StyleSelection } from 'sidebar/StyleSelection';
import { NETWORKS_DEFAULT_STYLE, NETWORK_STYLES } from 'config/networks/styles';
import { networksStyleState } from 'state/networks/networks-style';

export const NetworksSection: FC<{}> = () => {
  return (
    <SidebarPanel id="assets" title="Infrastructure">
      <Box p={2}>
        <NetworkControl />
      </Box>
      <Box p={2} bgcolor="#eee">
        <StyleSelection state={networksStyleState} options={NETWORK_STYLES} defaultValue={NETWORKS_DEFAULT_STYLE} />
      </Box>
    </SidebarPanel>
  );
};
