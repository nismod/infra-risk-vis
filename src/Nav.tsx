import { Menu } from '@mui/icons-material';
import {
  AppBar,
  Divider,
  Drawer,
  IconButton,
  List,
  ListItem,
  Link as MuiLink,
  Toolbar,
  Tooltip,
  styled,
} from '@mui/material';
import { Box } from '@mui/system';
import { FC, forwardRef, useCallback, useState } from 'react';
import { NavLink as RouterNavLink } from 'react-router-dom';

import { useIsMobile } from '../src/use-is-mobile';
import { withProps } from 'lib/react/with-props';
import { globalStyleVariables } from 'theme';

const Link = styled(MuiLink)({
  color: 'inherit',
  textDecoration: 'none',
});

const DrawerLink = styled(Link)({
  '&.active': {
    backgroundColor: '#eeeeee',
  },
});

const ToolbarLink = styled(Link)({
  padding: '0 3px 1px 3px',
  margin: '0 9px -10px 9px',
  borderBottom: '6px solid transparent',
  '&:hover,&:focus': {
    borderBottomColor: '#ffffffbc',
  },
  '&:active,&.active': {
    borderBottomColor: '#ffffff',
  },
});

const ToolbarNavLink = forwardRef<any, any>(({ ...others }, ref) => (
  <ToolbarLink variant="h6" component={RouterNavLink} ref={ref} {...others} />
));

const DrawerNavLink = forwardRef<any, any>(({ ...others }, ref) => (
  <DrawerLink component={RouterNavLink} ref={ref} {...others} />
));

const GrowingDivider = styled(Divider)({
  flexGrow: 1,
});

const drawerWidth = globalStyleVariables.mobileDrawerWidth;
const MobileDrawer = styled(Drawer)({
  width: drawerWidth,
  flexShrink: 0,
  [`& .MuiDrawer-paper`]: { width: drawerWidth, boxSizing: 'border-box' },
});

const NavTooltip = withProps(Tooltip, {
  enterDelay: 500,
  disableInteractive: true,
});

const MobileNavContent: FC<{ navItems: NavItemConfig[] }> = ({ navItems }) => {
  const [drawerOpen, setDrawerOpen] = useState(false);

  const closeDrawer = useCallback(() => {
    setDrawerOpen(false);
  }, []);

  return (
    <>
      <IconButton color="inherit" onClick={() => setDrawerOpen((open) => !open)} title="Menu">
        <Menu />
      </IconButton>

      <ToolbarNavLink to="/" onClick={closeDrawer}>
        J-SRAT
      </ToolbarNavLink>

      <GrowingDivider />

      <MobileDrawer open={drawerOpen} onClose={closeDrawer}>
        <Toolbar /> {/* Prevents app bar from concealing content*/}
        <List>
          <ListItem component={DrawerNavLink} to="/" exact onClick={closeDrawer}>
            Home
          </ListItem>
          {navItems.map(({ to, title, tooltip }) => (
            <NavTooltip key={to} title={tooltip} placement="right">
              <ListItem component={DrawerNavLink} to={to} onClick={closeDrawer}>
                {title}
              </ListItem>
            </NavTooltip>
          ))}
        </List>
      </MobileDrawer>
    </>
  );
};

const DesktopNavContent: FC<{ navItems: NavItemConfig[] }> = ({ navItems }) => (
  <>
    <ToolbarNavLink to="/">J-SRAT</ToolbarNavLink>

    {navItems.map(({ to, title, tooltip }) => (
      <NavTooltip title={tooltip} placement="bottom">
        <ToolbarNavLink key={to} to={to}>
          {title}
        </ToolbarNavLink>
      </NavTooltip>
    ))}
  </>
);

const topStripeHeight = 6;

export interface NavItemConfig {
  to: string;
  title: string;
  tooltip?: string;
}

export const Nav: FC<{ height: number; navItems: NavItemConfig[] }> = ({ height, navItems }) => {
  const isMobile = useIsMobile();

  return (
    <AppBar position="fixed" sx={{ color: 'white' }}>
      <Box height={topStripeHeight} width="100%" bgcolor="rgb(197,206,0)" />
      <Toolbar
        variant="dense"
        sx={{
          backgroundColor: 'rgb(0,126,133)',
          height: height - topStripeHeight,
        }}
      >
        {isMobile ? <MobileNavContent navItems={navItems} /> : <DesktopNavContent navItems={navItems} />}
      </Toolbar>
    </AppBar>
  );
};
