import React from 'react';
import { Alert, Box, Typography } from '@mui/material';
import { HazardsControl } from './controls/HazardsControl';
import { NetworkControl } from './controls/NetworkControl';
import { ViewModeToggle } from './ViewModeToggle';


export const SidebarContent = ({ view }) => {
  let component : React.ReactElement;
  switch(view){
    case 'exposure':
      component = (
        <>
          <NetworkControl />
          <Box mb={1}>
            <Typography variant="h6">View Mode</Typography>
            <ViewModeToggle />
          </Box>
          <HazardsControl />
        </>
      )
      break;
    case 'risk':
      component = (
        <Box mb={2}>
          <Typography variant="h6">Risk</Typography>
          <Alert severity="info">Coming soon.</Alert>
        </Box>);
      break;
    case 'adaptation':
      component = (
        <Box mb={2}>
          <Typography variant="h6">Adaptation</Typography>
          <Alert severity="info">Coming soon.</Alert>
        </Box>);
      break;
    case 'prioritization':
      component = (
        <Box mb={2}>
          <Typography variant="h6">Prioritization</Typography>
          <Alert severity="info">Coming soon.</Alert>
        </Box>);
      break;
    default:
      component = (<></>);
  }
  return component
};
