import React from 'react';
import { NavLink } from 'react-router-dom';
import { makeStyles } from '@material-ui/core/styles';
import AppBar from '@material-ui/core/AppBar';
import Toolbar from '@material-ui/core/Toolbar';
import Typography from '@material-ui/core/Typography';

const useStyles = makeStyles();

const Nav = () => {
  const classes = useStyles();
  return (
    <AppBar position="fixed">
      <Toolbar>
        <Typography variant="h6" className={classes.title}>
          <span className="nav-brand">
          Infrastructure Risk Assessment
          </span>
        </Typography>
        <NavLink exact className="nav-link" to='/'>
          <Typography variant="h6" className={classes.title}>
            About
          </Typography>
        </NavLink>
        <NavLink className="nav-link" to='/overview'>
          <Typography variant="h6" className={classes.title}>
            Infrastructure networks
          </Typography>
        </NavLink>
      </Toolbar>
    </AppBar>
  )
}

export default Nav
