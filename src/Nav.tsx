import React from 'react';
import { AppBar, Toolbar, Typography } from '@mui/material';
import { NavLink } from 'react-router-dom';

export const Nav = () => {
  return (
    <AppBar position="fixed">
      <Toolbar>
        <NavLink exact className="nav-link" to="/">
          <Typography variant="h6">J-SRAT</Typography>
        </NavLink>
        <NavLink className="nav-link" to="/exposure">
          <Typography variant="h6">Exposure</Typography>
        </NavLink>
        <NavLink className="nav-link" to="/risk">
          <Typography variant="h6">Risk</Typography>
        </NavLink>
        <NavLink className="nav-link" to="/adaptation">
          <Typography variant="h6">Adaptation</Typography>
        </NavLink>
        <NavLink className="nav-link" to="/prioritization">
          <Typography variant="h6">Prioritization</Typography>
        </NavLink>
        <NavLink className="nav-link" to="/data">
          <Typography variant="h6">Data</Typography>
        </NavLink>
      </Toolbar>
    </AppBar>
  );
};
