import React from 'react'
import { NavLink } from 'react-router-dom'

const Nav = () => (
  <nav className="navbar navbar-height navbar-expand navbar-dark">
    <NavLink className="navbar-brand" to="/">
      <img src="/logo.png" alt="OIA" />
    </NavLink>
    <ul className="navbar-nav mr-auto">
      <li className="nav-item">
        <NavLink className="nav-link" to='/overview'>
          Overview
        </NavLink>
      </li>
      <li className="nav-item">
        <NavLink className="nav-link" to='/flood'>
          Flood
        </NavLink>
      </li>
    </ul>
  </nav>
)

export default Nav
