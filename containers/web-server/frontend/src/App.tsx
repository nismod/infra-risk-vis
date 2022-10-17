import { CssBaseline, StyledEngineProvider, Toolbar } from '@mui/material';
import { ThemeProvider } from '@mui/material/styles';
import { Route, BrowserRouter as Router, Switch } from 'react-router-dom';
import { RecoilRoot } from 'recoil';

import { Nav } from './Nav';
import { DataPage } from './pages/DataPage';
import { IntroPage } from './pages/IntroPage';
import { MapPage } from './pages/MapPage';
import { theme } from './theme';

import './index.css';
import 'mapbox-gl/dist/mapbox-gl.css';

export const App = () => {
  return (
    <RecoilRoot>
      <StyledEngineProvider injectFirst>
        <ThemeProvider theme={theme}>
          <Router>
            <CssBaseline />
            <Nav />
            <Switch>
              <Route path="/" exact>
                <Toolbar /> {/* Prevents app bar from concealing content*/}
                <IntroPage />
              </Route>
              <Route path="/view/:view" render={({ match: { params } }) => <MapPage view={params.view} />} />
              <Route path="/data" exact>
                <Toolbar /> {/* Prevents app bar from concealing content*/}
                <DataPage />
              </Route>
            </Switch>
          </Router>
        </ThemeProvider>
      </StyledEngineProvider>
    </RecoilRoot>
  );
};
