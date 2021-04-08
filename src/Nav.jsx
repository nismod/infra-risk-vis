import React from 'react'
import { NavLink } from 'react-router-dom'

const Nav = () => (
  <nav className="navbar navbar-height navbar-expand navbar-light">
    <NavLink className="navbar-brand" to="/">
      <img src="/logo.png" alt="OIA" />
    </NavLink>
    <ul className="navbar-nav mr-auto">
      <li className="nav-item">
        <NavLink className="nav-link" to='/summary'>
          Summary
        </NavLink>
      </li>
      <li className="nav-item">
        <NavLink className="nav-link" to='/overview'>
          Overview
        </NavLink>
      </li>
      <li className="nav-item">
        <NavLink className="nav-link" to='/roads'>
          Roads
        </NavLink>
      </li>
      <li className="nav-item">
        <NavLink className="nav-link" to='/rail'>
          Rail
        </NavLink>
      </li>
      <li className="nav-item">
        <NavLink className="nav-link" to='/energy_network'>
          Electricity
        </NavLink>
      </li>
      <li className="nav-item">
        <NavLink className="nav-link" to='/flood'>
          Flood
        </NavLink>
      </li>
      <li className="nav-item">
        <NavLink className="nav-link" to='/risk'>
          Risk
        </NavLink>
      </li>
    </ul>
  </nav>
)

export default Nav
