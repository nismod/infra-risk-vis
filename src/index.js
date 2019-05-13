import React from 'react';
import { render } from 'react-dom'
import { Route, BrowserRouter as Router, Switch } from 'react-router-dom'

import Nav from './components/Nav'
import SimpleMap from './containers/SimpleMap'
import SelectMap from './containers/SelectMap'

import './index.css';
import 'bootstrap/dist/css/bootstrap.min.css'

render(
    <Router>
        <div>
            <Route path="/" component={Nav}/>
            <Switch>
                <Route path="/:name" render={props =>
                {
                    const {name} = props.match.params
                    switch(name) {
                        case 'overview':
                            return (
                                <SimpleMap
                                    key={'map_' + props.match.params.name}
                                    style={props.match.params.name}
                                />
                            )
                        case 'flood':
                            return (
                                <SelectMap
                                    key={'map_' + props.match.params.name}
                                    style={props.match.params.name}
                                />
                            )
                        default:
                            return (
                                <div>
                                    View not available
                                </div>
                            )
                    }
                }}/>
            </Switch>
        </div>
    </Router>,
    document.getElementById('root')
)
