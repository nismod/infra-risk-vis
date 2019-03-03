import React, { Component } from 'react'

import StaticMap from '../components/maps/StaticMap'
import TooltipMap from '../components/maps/TooltipMap'

const SimpleMap = ({ match }) => (
    <div>
        <TooltipMap style={"http://localhost:8080/styles/" + match.params.name + "/style.json"}/>
    </div>
)

export default SimpleMap
