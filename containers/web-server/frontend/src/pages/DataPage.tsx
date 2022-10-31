import { Alert, Paper, Table, TableBody, TableCell, TableContainer, TableHead, TableRow } from '@mui/material';
import React from 'react';

import ScrollToTop from '@/lib/react/ScrollToTop';

export const DataPage = () => (
  <article>
    <ScrollToTop />
    <h1>Data Sources and Access</h1>

    <h4>Disclaimer</h4>
    <p>
      This website was created for dissemination purposes. The information included here must not be used for the design
      of hazard-resistant structures or to support any important decision involving human life, capital or property. The
      values of risk or hazard in this map do not constitute an alternative nor do they replace building actions defined
      in national building codes. Readers seeking this information should consult national databases. The data presented
      here are the combination of results computed using multiple input models covering the majority of landmass. These
      models represent the best information currently publicly accessible, as far as the authors are aware. The
      presentation here results from an integration process that is solely the responsibility of the University of
      Oxford.
    </p>

    <Alert severity="info">TODO - Data comes from multiple sources</Alert>

    <Alert severity="info">TODO - Licensing</Alert>

    <h2>Hazard Data</h2>

    <Alert severity="info">TODO - Update Hazard Data</Alert>
    <TableContainer component={Paper}>
      <Table aria-label="simple table">
        <TableHead>
          <TableRow>
            <TableCell>Hazard type</TableCell>
            <TableCell>Data source</TableCell>
            <TableCell>Exceedance Probabilities (1/return periods in years)</TableCell>
            <TableCell>Intensities and spatial extents</TableCell>
            <TableCell>Climate scenario information</TableCell>
            <TableCell>Citation</TableCell>
            <TableCell>License</TableCell>
          </TableRow>
        </TableHead>
        <TableBody>
          <TableRow>
            <TableCell>Coastal flooding</TableCell>
            <TableCell>
              <a href="http://wri-projects.s3.amazonaws.com/AqueductFloodTool/download/v2/index.html">
                Aqueduct Floods Hazard Maps - Coastal
              </a>
            </TableCell>
            <TableCell>1/2, 1/5, 1/10, 1/25, 1/50, 1/100, 1/250, 1/500, 1/1000</TableCell>
            <TableCell>Flood depths in meters over 1km grid squares</TableCell>
            <TableCell>
              RCP Baseline, 4.5 &amp; 8.5 emission&nbsp;scenarios.
              <br />
              Current + future maps in 2030, 2050 and 2080
            </TableCell>
          </TableRow>
          <TableRow>
            <TableCell>Riverine Flooding</TableCell>
            <TableCell>
              <a href="http://wri-projects.s3.amazonaws.com/AqueductFloodTool/download/v2/index.html">
                Aqueduct Floods Hazard Maps - Riverine
              </a>
            </TableCell>
            <TableCell>1/2, 1/5, 1/10, 1/25, 1/50, 1/100, 1/250, 1/500, 1/1000</TableCell>
            <TableCell>Flood depths in meters over 1km grid squares </TableCell>
            <TableCell>
              RCP Baseline, 4.5 &amp; 8.5 emission&nbsp;scenarios.
              <br />
              Current + future maps in 2030, 2050 and 2080
            </TableCell>
          </TableRow>
          <TableRow>
            <TableCell>Extreme Heat (Exposure)</TableCell>
            <TableCell>
              <a href="https://data.isimip.org/datasets/4f79e7aa-7def-4665-854f-93ff033bec37/">ISIMP Extreme Heat</a>
            </TableCell>
            <TableCell>TODO</TableCell>
            <TableCell>TODO</TableCell>
            <TableCell>
              RCP Baseline, 2.6 &amp; 6.0 emission&nbsp;scenarios.
              <br />
              Current + future maps in 2030, 2050 and 2080
            </TableCell>
          </TableRow>
          <TableRow>
            <TableCell>STORM (Cyclones)</TableCell>
            <TableCell>
              <a href="https://data.4tu.nl/articles/dataset/STORM_tropical_cyclone_wind_speed_return_periods/12705164/3">
                STORM Tropical Cyclone
              </a>
            </TableCell>
            <TableCell>
              10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000,
              5000, 6000, 7000, 8000, 9000, 10000, 11000
            </TableCell>
            <TableCell>Wind Speeds at Fixed Return Periods</TableCell>
            <TableCell>
              RCP Baseline
              <br />
              Current + future maps in 2030, 2050 and 2080
            </TableCell>
          </TableRow>
          <TableRow>
            <TableCell>Seismic Risk</TableCell>
            <TableCell>TODO</TableCell>
            <TableCell></TableCell>
            <TableCell></TableCell>
            <TableCell></TableCell>
            <TableCell>
              Pagani M, Garcia-Pelaez J, Gee R, Johnson K, Silva V, Simionato M, Styron R, Vigano D, Danciu L, Monelli
              D, Poggi V, Weatherill G. (2019). The 2018 version of the Global Earthquake Model: Hazard component.
              Earthquake Spectra, 36(1), https://doi.org/10.1177/87552930209318.
            </TableCell>
            <TableCell>
              The visualisation of seismic hazard data is licensed under the terms of the{' '}
              <a href="https://creativecommons.org/licenses/by-nc-sa/4.0/">
                Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License (CC BY-NC-SA)
              </a>
              . The underlying data are confidential, contact GEM for more information.
            </TableCell>
          </TableRow>
        </TableBody>
      </Table>
    </TableContainer>

    <h2>Infrastructure Network Data</h2>

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
