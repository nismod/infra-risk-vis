import { AppBar, Toolbar, Typography } from '@mui/material';
import React from 'react';
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
          <Typography variant="h6">Global</Typography>
        </NavLink>
        <NavLink className="nav-link" to="/view/hazard">
          <Typography variant="h6">Hazard</Typography>
        </NavLink>
        <NavLink className="nav-link" to="/view/exposure">
          <Typography variant="h6">Exposure</Typography>
        </NavLink>
        <NavLink className="nav-link" to="/view/vulnerability">
          <Typography variant="h6">Vulnerability</Typography>
        </NavLink>
        <NavLink className="nav-link" to="/view/risk">
          <Typography variant="h6">Risk</Typography>
        </NavLink>
        <NavLink className="nav-link" to="/data">
          <Typography variant="h6">Data</Typography>
        </NavLink>
      </Toolbar>
    </AppBar>
  );
};
