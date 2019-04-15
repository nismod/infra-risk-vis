import React from 'react'

import StaticMap from '../components/maps/StaticMap'

const SimpleMap = ({ style }) => (
    <div className="d-flex align-items-stretch map-height">
        <div className="col-sm-12">
            <StaticMap 
                style={"http://localhost:8080/styles/" + style + "/style.json"}
                toggleableLayerIds={['road_secondary_tertiary', 'road_trunk_primary', 'road_major_motorway', 'aeroway', 'waterway']}/>
        </div>
    </div>
)

export default SimpleMap
