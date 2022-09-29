import React from 'react';
import { AppBar, Toolbar, Typography } from '@mui/material';
import { NavLink } from 'react-router-dom';

export const Nav = () => {
  return (
    <AppBar position="fixed">
      <Toolbar
        sx={{
          background:
            'linear-gradient(180deg, rgba(0,128,0,1) 0%, rgba(0,128,0,1) 5%, rgba(9,9,9,1) 5%, rgba(9,9,9,1) 100%);',
          '& a.nav-link': {
            '&:hover,&:focus,&:active,&.active': {
              borderBottom: '1px solid #ffd700',
            },
          },
        }}
      >
        <NavLink exact className="nav-link" to="/">
          <Typography variant="h6">G-SRAT</Typography>
        </NavLink>
        <NavLink className="nav-link" to="/exposure">
          <Typography variant="h6">Exposure</Typography>
        </NavLink>
        {/* <NavLink className="nav-link" to="/risk">
          <Typography variant="h6">Risk</Typography>
        </NavLink>
        <NavLink className="nav-link" to="/adaptation">
          <Typography variant="h6">Adaptation</Typography>
        </NavLink>
        <NavLink className="nav-link" to="/nature-based-solutions">
          <Typography variant="h6">Nature-based Solutions</Typography>
        </NavLink> */}
        <NavLink className="nav-link" to="/data">
          <Typography variant="h6">Data</Typography>
        </NavLink>
      </Toolbar>
    </AppBar>
  );
};
