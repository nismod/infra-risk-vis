import React from 'react'

import StaticMap from '../components/maps/StaticMap'

const SimpleMap = ({ style }) => (
    <div className="map-height">
        <StaticMap
            style={"http://localhost:8080/styles/" + style + "/style.json"}
            toggleableLayerIds={['road_secondary_tertiary', 'road_trunk_primary', 'road_major_motorway', 'aeroway', 'waterway']}
            clickableLayerAttributes={{
                'road_secondary_tertiary': {
                    '_header': 'Rural Road',
                    'gml_id': 'Name',
                    'u_jurisdic': 'Classification'
                },
                'road_trunk_primary': {
                    '_header': 'Provincial Road',
                    'name': 'Name',
                    'jurisdicti': 'Classification'
                },
                'road_major_motorway': {
                    '_header': 'Motorway',
                    'name': 'Name',
                    'jurisdicti': 'Classification'
                },
                'waterway': {},
                'aeroway': {
                    '_header': 'Flight path',
                    'from_iata': 'From',
                    'to_iata': 'To'
                }
            }}/>
    </div>
)

export default SimpleMap
