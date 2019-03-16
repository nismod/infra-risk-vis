import React, {Component} from 'react'
import HighlightMap from '../components/maps/HighlightMap'

class AttributesMap extends Component {
    constructor(props) {
      super(props);
      this.onHighlight = this.onHighlight.bind(this);
      this.state = {
          highFeatures: []
      }
    }

    onHighlight(features) {
      this.setState({
        highFeatures: features
      });
    }

    render() {
      return (
        <div className="d-flex align-items-stretch map-height">
            <div className="col-sm-8">
                <HighlightMap 
                    style={"http://localhost:8080/styles/" + this.props.style + "/style.json"}
                    enabledFeatures={['road_edges_national', 'road_edges_provincial', 'road_edges_rural', 'water_edges', 'water_nodes']}
                    onHighlight={this.onHighlight}/>
            </div>

            <div className="col-sm-4">
                {this.state.highFeatures}
            </div>
        </div>
      )
    }
  }

export default AttributesMap
