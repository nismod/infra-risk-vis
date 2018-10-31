import React from 'react';
import { render } from 'react-dom'
import { Route, BrowserRouter as Router, Switch } from 'react-router-dom'

import Nav from './components/Nav'
import SimpleMap from './containers/SimpleMap'
import * as serviceWorker from './serviceWorker';

import './index.css';
import 'bootstrap/dist/css/bootstrap.min.css'

render(
    <Router>
        <div>
            <Route path="/" component={Nav}/>
            <main role="main" className="container">
                <Switch>    
                    <Route path="/:name" render={props =>
                        <SimpleMap 
                            key={props.match.params.name}
                            match={props.match}
                        />}/>
                </Switch>
            </main>
        </div>
    </Router>,
    document.getElementById('root')
)
// If you want your app to work offline and load faster, you can change
// unregister() to register() below. Note this comes with some pitfalls.
// Learn more about service workers: http://bit.ly/CRA-PWA
serviceWorker.unregister();
