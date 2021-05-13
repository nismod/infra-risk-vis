import React from 'react'
import { NavLink } from 'react-router-dom'

const Nav = () => (
  <nav className="navbar navbar-height navbar-expand navbar-light">
    <ul className="navbar-nav mr-auto">
      <li className="nav-item">
        <NavLink exact className="nav-link" to='/'>
          About
        </NavLink>
      </li>
      <li className="nav-item">
        <NavLink className="nav-link" to='/overview'>
          Infrastructure networks
        </NavLink>
      </li>
    </ul>
  </nav>
)

export default Nav
