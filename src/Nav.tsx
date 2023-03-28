import React from 'react';
import { AppBar, Toolbar, Tooltip, Typography } from '@mui/material';
import { NavLink } from 'react-router-dom';

export const Nav = () => {
  return (
    <AppBar position="fixed">
      <Toolbar
        sx={{
          background:
            'linear-gradient(180deg, rgba(197,206,0,1) 0%, rgba(197,206,0,1) 10%, rgba(0,126,133,1) 10%, rgba(0,126,133,1) 100%);',
          '& a.nav-link': {
            '&:hover,&:focus,&:active,&.active': {
              borderBottomColor: '#ffffff',
            },
          },
        }}
      >
        <NavLink exact className="nav-link" to="/">
          <Typography variant="h6">J-SRAT</Typography>
        </NavLink>
        <NavTooltip title="Infrastructure assets and natural hazards">
          <NavLink className="nav-link" to="/exposure">
            <Typography variant="h6">Exposure</Typography>
          </NavLink>
        </NavTooltip>
        <NavTooltip title="Risk of hazard-related damages to assets">
          <NavLink className="nav-link" to="/risk">
            <Typography variant="h6">Risk</Typography>
          </NavLink>
        </NavTooltip>
        <NavTooltip title="Adaptation options to decrease hazard-related risks">
          <NavLink className="nav-link" to="/adaptation">
            <Typography variant="h6">Adaptation</Typography>
          </NavLink>
        </NavTooltip>
        <NavTooltip title="Analysis of nature-based solutions potential">
          <NavLink className="nav-link" to="/nature-based-solutions">
            <Typography variant="h6">Nature-based Solutions</Typography>
          </NavLink>
        </NavTooltip>
        <NavTooltip title="More information about datasets in the tool">
          <NavLink className="nav-link" to="/data">
            <Typography variant="h6">About</Typography>
          </NavLink>
        </NavTooltip>
      </Toolbar>
    </AppBar>
  );
};

function NavTooltip({ ...props }: any) {
  return <Tooltip enterDelay={500} {...props} />;
}
