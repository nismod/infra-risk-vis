import React from 'react'
import { NavLink } from 'react-router-dom'

const Nav = () => (
  <nav className="navbar navbar-height navbar-expand navbar-dark">
    <NavLink className="navbar-brand" to="/">
      <img src="/logo.png" alt="OIA" />
    </NavLink>
    <ul className="navbar-nav mr-auto">
      <li className="nav-item">
        <a className="nav-link" href='/overview'>
          Overview
        </a>
      </li>
      <li className="nav-item">
        <a className="nav-link" href='/roads'>
          Roads
        </a>
      </li>
      <li className="nav-item">
        <a className="nav-link" href='/energy_network'>
          Energy Network
        </a>
      </li>
      <li className="nav-item">
        <a className="nav-link" href='/flood'>
          Flood
        </a>
      </li>
      <li className="nav-item">
        <a className="nav-link" href='/impact'>
          Impact
        </a>
      </li>
      <li className="nav-item">
        <a className="nav-link" href='/risk'>
          Risk
        </a>
      </li>
      <li className="nav-item">
        <a className="nav-link" href='/adaptation'>
          Assets with Adaptation BCR &gt; 1
        </a>
      </li>
    </ul>
  </nav>
)

export default Nav
