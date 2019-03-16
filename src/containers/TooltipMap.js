import React from 'react'

// import StaticMap from '../components/maps/StaticMap'
import TooltipMap from '../components/maps/TooltipMap'

const SimpleMap = ({ style }) => (
    <div className="d-flex align-items-stretch map-height">
        <div className="col-sm-12">
            <TooltipMap style={"http://localhost:8080/styles/" + style + "/style.json"}/>
        </div>
    </div>
)

export default SimpleMap
