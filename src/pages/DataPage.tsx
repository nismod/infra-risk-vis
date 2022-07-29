import React from 'react';
import { Alert, Button, Paper, Table, TableBody, TableCell, TableContainer, TableHead, TableRow, Typography } from '@mui/material';
import ScrollToTop from 'lib/hooks/scroll-to-top';

export const DataPage = () => (
  <article>
    <ScrollToTop />
    <p>
      <Alert
        severity="success"
        action={
          <Button color="inherit" size="small">
            <a href="https://github.com/nismod/infra-risk-vis/issues">
              REPORT
            </a>
          </Button>
        }
      >
        The tool has recently been released. Please tell us if anything is not
        working as it should and suggest potential improvements.
      </Alert>
    </p>

    <p>The modelling and analysis presented here aim to support climate
      adaptation decision-making by identifying spatial criticalities and risks
      under current and future climate scenarios.</p>

    <p>The following table summarises the data and model results presented.</p>

    <TableContainer component={Paper}>
      <Table aria-label="simple table">
        <TableHead>
          <TableRow>
            <TableCell>Infrastructure Sector</TableCell>
            <TableCell>Assets</TableCell>
            <TableCell>Exposure</TableCell>
          </TableRow>
        </TableHead>
        <TableBody>
          <TableRow>
            <TableCell>Transport</TableCell>
            <TableCell>Road links</TableCell>
            <TableCell>Lengths exposed to coastal or river flooding</TableCell>
          </TableRow>
          <TableRow>
            <TableCell>Energy</TableCell>
            <TableCell>Electricity generation (power stations) and transmission lines</TableCell>
            <TableCell>Lengths exposed to high wind speeds from tropical cyclones</TableCell>
          </TableRow>
        </TableBody>
      </Table>
    </TableContainer>

    <p>For more details on the infrastructure and hazard data used, see below</p>

    <p>The primary output metrics from the analysis are:</p>

    <ul>
      <li>Exposure of infrastructure assets to flooding and cyclone winds</li>
      <li>Expected Annual Damages (EAD) (direct physical risks) estimated as the
        area under the direct damage vs exceedance probability curve </li>
    </ul>


    <h2>Open-source code</h2>

    <p>This tool to visualize the model outputs is developed and documented
      here:</p>

    <ul>
      <li>
        <a href="https://github.com/nismod/infra-risk-vis" target="blank">
          github.com/nismod/infra-risk-vis
        </a>
      </li>
    </ul>

    <p>The analytics for the Caribbean are produced using the code and models here:</p>

    <ul>
      <li>
        <a href="https://github.com/nismod/open-gira" target="blank">
          github.com/nismod/open-gira
        </a>
      </li>
    </ul>

    <h1>Data Sources and Access</h1>

    <p>Data comes from multiple open data sources.</p>

    <h2>Hazard Data</h2>

    <TableContainer component={Paper}>
      <Table aria-label="simple table">
        <TableHead>
          <TableRow>
            <TableCell>Hazard type</TableCell>
            <TableCell>Data source</TableCell>
            <TableCell>Exceedance Probabilities (1/return periods in years)</TableCell>
            <TableCell>Intensities and spatial extents</TableCell>
            <TableCell>Climate scenario information</TableCell>
          </TableRow>
        </TableHead>
        <TableBody>
          <TableRow>
            <TableCell>River and coastal flooding</TableCell>
            <TableCell>WRI Aqueduct</TableCell>
            <TableCell>1/20, 1/50, 1/100, 1/200, 1/500, and 1/1500</TableCell>
            <TableCell>Flood depths in meters over 100m grid squares</TableCell>
            <TableCell>RCP 2.6, 4.5 &amp; 8.5 emission&nbsp;scenarios.<br/>Current + future maps in 2050 and 2080</TableCell>
          </TableRow>
          <TableRow>
            <TableCell>Tropical cyclones (winds)</TableCell>
            <TableCell><a href="https://data.4tu.nl/articles/dataset/STORM_IBTrACS_present_climate_synthetic_tropical_cyclone_tracks/12706085/2">STORM IBTrACS model</a></TableCell>
            <TableCell>Annual exceedance probabilities from 1/10 to 1/1000</TableCell>
            <TableCell>10 minute sustained maximum wind speeds in m/s at 10km grid squares</TableCell>
            <TableCell>RCP 8.5 emission scenarios.<br/>Current + future maps in 2050</TableCell>
          </TableRow>
        </TableBody>
      </Table>
    </TableContainer>


    <h2>Infrastructure Network Data</h2>

    <TableContainer component={Paper}>
      <Table aria-label="simple table">
        <TableHead>
          <TableRow>
            <TableCell>Sector</TableCell>
            <TableCell>Sub-sector</TableCell>
            <TableCell>Asset highlights</TableCell>
            <TableCell>Data sources</TableCell>
          </TableRow>
        </TableHead>
        <TableBody>
          <TableRow>
            <TableCell rowSpan={2}>Energy</TableCell>
            <TableCell>Generation</TableCell>
            <TableCell>power plants</TableCell>
            <TableCell>Global Powerplants Database</TableCell>
          </TableRow>
          <TableRow>
            <TableCell>Transmission lines</TableCell>
            <TableCell>km lines</TableCell>
            <TableCell>Gridfinder</TableCell>
          </TableRow>
          <TableRow>
            <TableCell rowSpan={2}>Transport</TableCell>
            <TableCell>Railways</TableCell>
            <TableCell>stations, km tracks</TableCell>
          </TableRow>
          <TableRow>
            <TableCell>Roads</TableCell>
            <TableCell>km roads</TableCell>
          </TableRow>
          <TableRow>
            <TableCell>Buildings</TableCell>
            <TableCell>buildings</TableCell>
            <TableCell>OpenStreetMap</TableCell>
          </TableRow>
        </TableBody>
      </Table>
    </TableContainer>

    <h2>Contextual Map Data</h2>

    <p>Background map data is &copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors, style &copy; <a href="https://carto.com/attributions">CARTO</a>.</p>

    <p>Satellite imagery background is derived from <a href="https://s2maps.eu">Sentinel-2 cloudless - https://s2maps.eu</a> by <a href="https://eox.at">EOX IT Services GmbH</a> (Contains modified Copernicus Sentinel data 2020).</p>
  </article>
);
