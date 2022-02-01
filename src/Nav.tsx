import React from 'react';
import { AppBar, Toolbar, Typography } from '@mui/material';
import { NavLink } from 'react-router-dom';

export const Nav = () => {
  return (
    <AppBar position="fixed">
      <Toolbar>
        <Typography variant="h6">
          <span className="nav-brand">Infrastructure Risk Assessment</span>
        </Typography>
        <NavLink exact className="nav-link" to="/">
          <Typography variant="h6">About</Typography>
        </NavLink>
        <NavLink className="nav-link" to="/overview">
          <Typography variant="h6">Infrastructure networks</Typography>
        </NavLink>
        {/* <NavLink className="nav-link" to="/hazards">
          <Typography variant="h6">Hazards</Typography>
        </NavLink> */}
      </Toolbar>
    </AppBar>
  );
};
