import React, { Component } from 'react'
import { NavLink } from 'react-router-dom'

class Nav extends Component {

    render() {
        return (
            <nav className="navbar navbar-height navbar-expand-md navbar-dark">
                <a className="navbar-brand" href="/"><img src="/logo.png" alt="OIA" /></a>
                <div className="collapse navbar-collapse" id="navbarCollapse">
                    <ul className="navbar-nav mr-auto">
                        <li className="nav-item">
                            <NavLink
                                className="nav-link"
                                to={'/overview'}>
                                Overview
                            </NavLink>
                        </li>
                        <li className="nav-item">
                            <NavLink
                                className="nav-link"
                                to={'/flood'}>
                                Flood
                            </NavLink>
                        </li>
                    </ul>
                </div>
            </nav>
        )
    }
}

export default Nav
