import React from 'react'

import TooltipMap from '../components/maps/TooltipMap'

const SelectMap = ({ style }) => (
    <div className="map-height">
        <TooltipMap
            style={"http://localhost:8080/styles/" + style + "/style.json"}
            tooltipLayerSources={['flood']}/>
    </div>
)

export default SelectMap
