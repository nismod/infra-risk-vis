import React from 'react';
import { Alert, Button, Paper, Table, TableBody, TableCell, TableContainer, TableHead, TableRow } from '@mui/material';
import ScrollToTop from 'lib/hooks/scroll-to-top';

export const DataPage = () => (
  <article>
    <ScrollToTop />
    <p>
      <Alert
        severity="success"
        action={
          <Button color="inherit" size="small">
            <a href="https://github.com/nismod/infra-risk-vis/issues">REPORT</a>
          </Button>
        }
      >
        The tool has recently been released. Please tell us if anything is not working as it should and suggest
        potential improvements.
      </Alert>
    </p>

    <p>
      The modelling and analysis presented here aim to support climate adaptation decision-making by identifying spatial
      criticalities and risks under current and future climate scenarios.
    </p>

    <p>The following table summarises the data and model results presented.</p>

    <TableContainer component={Paper} sx={{ my: 2 }}>
      <Table aria-label="simple table">
        <TableHead>
          <TableRow>
            <TableCell>Infrastructure Sector</TableCell>
            <TableCell>Assets</TableCell>
            <TableCell>Expected Annual Damages (EAD)</TableCell>
            <TableCell>Expected Annual Economic Losses (EAEL)</TableCell>
          </TableRow>
        </TableHead>
        <TableBody>
          <TableRow>
            <TableCell>Road</TableCell>
            <TableCell>Trunk, motorway, primary, secondary and tertiary roads</TableCell>
            <TableCell>Cost of rehabilitation/reinstating damaged assets</TableCell>
            <TableCell>Rerouting costs + wider effects of service disruption</TableCell>
          </TableRow>
          <TableRow>
            <TableCell>Rail</TableCell>
            <TableCell>Rail lines and stations</TableCell>
            <TableCell>Cost of rehabilitation/reinstating damaged assets</TableCell>
            <TableCell>Wider effects of service disruption</TableCell>
          </TableRow>
          <TableRow>
            <TableCell>Ports and airports</TableCell>
            <TableCell>Major sea and inland ports, international airports</TableCell>
            <TableCell>Not assessed for damages</TableCell>
            <TableCell>
              Included in the network as sources/sinks for transport flow mapping, but not assessed for effects of
              service disruption due to flooding.
            </TableCell>
          </TableRow>
        </TableBody>
      </Table>
    </TableContainer>

    <p>For more details on the infrastructure and hazard data used, see below</p>

    <p>The primary output metrics from the analysis are:</p>

    <ul>
      <li>
        Expected Annual Damages (EAD) (direct physical risks) estimated as the area under the direct damage vs
        exceedance probability curve{' '}
      </li>
      <li>
        Expected Annual Economic Losses (EAEL) (indirect economic risks) estimated as the area under the economic loss
        vs exceedance probability curve{' '}
      </li>
    </ul>

    <h2>Open-source code</h2>

    <p>This tool to visualize the model outputs is developed and documented here:</p>

    <ul>
      <li>
        <a href="https://github.com/nismod/infra-risk-vis" target="blank">
          github.com/nismod/infra-risk-vis
        </a>
      </li>
    </ul>

    <p>The analytics for Kenya, Tanzania, Uganda and Zambia are produced using the code and models here:</p>

    <ul>
      <li>
        <a href="https://github.com/nismod/east-africa-transport" target="blank">
          github.com/nismod/east-africa-transport
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
            <TableCell>
              <a href="https://www.wri.org/data/aqueduct-floods-hazard-maps">WRI Aqueduct Floods Hazard Maps</a>
            </TableCell>
            <TableCell>1/20, 1/50, 1/100, 1/200, 1/500, and 1/1500</TableCell>
            <TableCell>Flood depths in meters on a ~1km grid</TableCell>
            <TableCell>
              Current climate and future RCP&nbsp;4.5 and 8.5 emission scenarios in 2030, 2050 and 2080
            </TableCell>
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
            <TableCell>Transport</TableCell>
            <TableCell>Roads</TableCell>
            <TableCell>km roads</TableCell>
            <TableCell>
              <a href="https://www.openstreetmap.org/#map=6/-6.599/32.278" target="_blank" rel="noopener noreferrer">
                OpenStreetMap
              </a>
            </TableCell>
          </TableRow>
          <TableRow>
            <TableCell>Transport</TableCell>
            <TableCell>Rail</TableCell>
            <TableCell>km rail lines, stations</TableCell>
            <TableCell>
              <a href="https://www.openstreetmap.org/#map=6/-6.599/32.278" target="_blank" rel="noopener noreferrer">
                OpenStreetMap
              </a>
            </TableCell>
          </TableRow>
        </TableBody>
      </Table>
    </TableContainer>

    <h2>Contextual Map Data</h2>

    <p>
      Background map data is &copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors,
      style &copy; <a href="https://carto.com/attributions">CARTO</a>.
    </p>

    <p>
      Satellite imagery background is derived from{' '}
      <a href="https://s2maps.eu">Sentinel-2 cloudless - https://s2maps.eu</a> by{' '}
      <a href="https://eox.at">EOX IT Services GmbH</a> (Contains modified Copernicus Sentinel data 2020).
    </p>
  </article>
);
