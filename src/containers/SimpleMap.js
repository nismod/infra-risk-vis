import React, { Component } from 'react'

import StaticMap from '../components/Maps'

const SimpleMap = ({ match }) => (
    <div>
        <StaticMap style={"http://localhost:8080/styles/" + match.params.name + "/style.json"}/>
    </div>
)

export default SimpleMap
