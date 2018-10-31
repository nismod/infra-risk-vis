import React, { Component } from 'react'
import { NavLink, Route } from 'react-router-dom'

class Nav extends Component {

    render() {
        return (
            <nav className="navbar navbar-expand-md navbar-dark fixed-top bg-dark">
                <a className="navbar-brand" href="#">OI Analytics</a>
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
                                to={'/road'}>
                                Roads
                            </NavLink>
                        </li>
                        <li className="nav-item">
                            <NavLink
                                className="nav-link"
                                to={'/water'}>
                                Water
                            </NavLink>
                        </li>
                        <li className="nav-item">
                            <NavLink
                                className="nav-link"
                                to={'/air'}>
                                Air
                            </NavLink>
                        </li>
                    </ul>
                </div>
            </nav>
        )
    }
}

export default Nav
