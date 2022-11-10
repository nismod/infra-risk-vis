import React, { useState } from 'react';
import { AppBar, Divider, Drawer, IconButton, MenuList, MenuItem, Toolbar, Typography, useMediaQuery } from '@mui/material';
import { Menu } from "@mui/icons-material";
import { NavLink } from 'react-router-dom';

import { globalStyleVariables } from '@/theme';

export const Nav = () => {
  const isMobile = useMediaQuery((theme: any) => theme.breakpoints.down('md'));
  const [openDrawer, setOpenDrawer] = useState(false);


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
          minHeight: globalStyleVariables.navbarHeight
        }}
      >
        {
          isMobile ? (
            <>
              <IconButton
                onClick={() => setOpenDrawer(!openDrawer)}
                sx={{mt:1, pl:0}}
                title="Menu"
                >
                <Menu htmlColor='white' />
              </IconButton>
              <NavLink exact className="nav-link" to="/" onClick={() => setOpenDrawer(false)}>
                <Typography variant="h6">G-SRAT</Typography>
              </NavLink>
              <Divider sx={{flexGrow: 1}} />
              <a className="nav-link" href="http://www.globalresilienceindex.org/" target="_blank" rel="noopener noreferrer">
                <img
                  height="35"
                  src="/logo-grii-white.png"
                  alt="GRII"
                  />
              </a>

              <Drawer
                open={openDrawer}
                onClose={() => setOpenDrawer(false)}
                sx={{width: 240}}
              >
                <Toolbar sx={{width: 240}} /> {/* Prevents app bar from concealing content*/}
                <MenuList>
                  <MenuItem sx={{p:0}} onClick={() => setOpenDrawer(false)}>
                    <NavLink exact className="nav-link in-drawer" to="/">
                      Home
                    </NavLink>
                  </MenuItem>
                  <MenuItem sx={{p:0}} onClick={() => setOpenDrawer(false)}>
                    <NavLink className="nav-link in-drawer" to="/view/hazard">
                      Hazard
                    </NavLink>
                  </MenuItem>

                  <MenuItem sx={{p:0}} onClick={() => setOpenDrawer(false)}>
                    <NavLink className="nav-link in-drawer" to="/view/exposure">
                      Exposure
                    </NavLink>
                  </MenuItem>

                  <MenuItem sx={{p:0}} onClick={() => setOpenDrawer(false)}>
                    <NavLink className="nav-link in-drawer" to="/view/vulnerability">
                      Vulnerability
                    </NavLink>
                  </MenuItem>

                  <MenuItem sx={{p:0}} onClick={() => setOpenDrawer(false)}>
                    <NavLink className="nav-link in-drawer" to="/view/risk">
                      Risk
                    </NavLink>
                  </MenuItem>

                  <MenuItem sx={{p:0}} onClick={() => setOpenDrawer(false)}>
                    <NavLink className="nav-link in-drawer" to="/data">
                      About
                    </NavLink>
                  </MenuItem>
                </MenuList>
              </Drawer>
            </>
          ) : (
            <>
              <NavLink exact className="nav-link" to="/">
                <Typography variant="h6">G-SRAT</Typography>
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
                <Typography variant="h6">About</Typography>
              </NavLink>
              <Divider sx={{flexGrow: 1}} />
              <a className="nav-link" href="http://www.globalresilienceindex.org/" target="_blank" rel="noopener noreferrer">
                <img
                  height="35"
                  src="/logo-grii-white.png"
                  alt="GRII"
                />
              </a>
            </>
          )
        }
      </Toolbar>
    </AppBar>
  );
};
