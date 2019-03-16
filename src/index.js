import React from 'react';
import { render } from 'react-dom'
import { Route, BrowserRouter as Router, Switch } from 'react-router-dom'

import Nav from './components/Nav'
import SimpleMap from './containers/SimpleMap'
import TooltipMap from './containers/TooltipMap'
import AttributesMap from './containers/AttributesMap'
import * as serviceWorker from './serviceWorker';

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
                        case 'road':
                            return (
                                <AttributesMap 
                                    key={'map_' + props.match.params.name}
                                    style={props.match.params.name}
                                />
                            )
                        case 'water':
                            return (
                                <AttributesMap 
                                    key={'map_' + props.match.params.name}
                                    style={props.match.params.name}
                                />
                            )
                        case 'air':
                            return (
                                <TooltipMap 
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
// If you want your app to work offline and load faster, you can change
// unregister() to register() below. Note this comes with some pitfalls.
// Learn more about service workers: http://bit.ly/CRA-PWA
serviceWorker.unregister();
