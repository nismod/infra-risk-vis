import React from 'react'

import TooltipMap from '../components/maps/TooltipMap'

const SelectMap = ({ style }) => (
    <div className="d-flex align-items-stretch map-height">
        <div className="col-sm-12">
            <TooltipMap
                style={"http://localhost:8080/styles/" + style + "/style.json"}
                tooltipLayerSources={['flood']}/>
        </div>
    </div>
)

export default SelectMap
