import React from 'react';
import { Alert, Paper, Table, TableBody, TableCell, TableContainer, TableHead, TableRow } from '@mui/material';
import ScrollToTop from 'lib/hooks/scroll-to-top';

export const DataPage = () => (
  <article>
    <ScrollToTop />
    <h1>Data Sources and Access</h1>

    <p>TODO - Data comes from multiple sources</p>

    <Alert severity="info">TODO - Licensing</Alert>

    <h2>TODO - Update Hazard Data</h2>

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
            <TableCell>Coastal flooding</TableCell>
            <TableCell>
              <a href="http://wri-projects.s3.amazonaws.com/AqueductFloodTool/download/v2/index.html">Aqueduct Floods Hazard Maps - Coastal</a></TableCell>
            <TableCell>1/2, 1/5, 1/10, 1/25, 1/50, 1/100, 1/250, 1/500, 1/1000</TableCell>
            <TableCell>Flood depths in meters over 1km grid squares</TableCell>
            <TableCell>RCP Baseline, 4.5 &amp; 8.5 emission&nbsp;scenarios.<br />Current + future maps in 2030, 2050 and 2080</TableCell>
          </TableRow>
          <TableRow>
            <TableCell>Riverine Flooding</TableCell>
            <TableCell><a href="http://wri-projects.s3.amazonaws.com/AqueductFloodTool/download/v2/index.html">Aqueduct Floods Hazard Maps - Riverine</a></TableCell>
            <TableCell>1/2, 1/5, 1/10, 1/25, 1/50, 1/100, 1/250, 1/500, 1/1000</TableCell>
            <TableCell>Flood depths in meters over 1km grid squares </TableCell>
            <TableCell>RCP Baseline, 4.5 &amp; 8.5 emission&nbsp;scenarios.<br />Current + future maps in 2030, 2050 and 2080</TableCell>
          </TableRow>
        </TableBody>
      </Table>
    </TableContainer>


    <h2>Infrastructure Network Data</h2>

    <h2>Contextual Map Data</h2>

    <p>Background map data is &copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors, style &copy; <a href="https://carto.com/attributions">CARTO</a>.</p>

    <p>Satellite imagery background is derived from <a href="https://s2maps.eu">Sentinel-2 cloudless - https://s2maps.eu</a> by <a href="https://eox.at">EOX IT Services GmbH</a> (Contains modified Copernicus Sentinel data 2020).</p>
  </article>
);