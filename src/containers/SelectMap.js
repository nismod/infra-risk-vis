import React from 'react'

import StaticMap from '../components/maps/StaticMap'

const SelectMap = ({ style }) => (
    <div className="d-flex align-items-stretch map-height">
        <div className="col-sm-12">
            <StaticMap style={"http://localhost:8080/styles/" + style + "/style.json"}/>
        </div>
    </div>
)

export default SelectMap
