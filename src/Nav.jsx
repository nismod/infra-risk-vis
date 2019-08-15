import React from 'react'
import { NavLink } from 'react-router-dom'

const Nav = () => (
  <nav className="navbar navbar-height navbar-expand navbar-dark">
    <NavLink className="navbar-brand" to="/">
      <img src="/logo.png" alt="OIA" />
    </NavLink>
    <ul className="navbar-nav mr-auto">
      <li className="nav-item">
        <a className="nav-link" href='/roads'>
          Roads
        </a>
      </li>
      <li className="nav-item">
        <a className="nav-link" href='/rail'>
          Railways
        </a>
      </li>
      <li className="nav-item">
        <a className="nav-link" href='/airwater'>
          Airports and Waterway ports
        </a>
      </li>
      <li className="nav-item">
        <a className="nav-link" href='/flood'>
          Flood
        </a>
      </li>
    </ul>
  </nav>
)

export default Nav
