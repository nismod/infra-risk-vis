import React from 'react';
import { Alert, Button, Paper, Table, TableBody, TableCell, TableContainer, TableHead, TableRow } from '@mui/material';
import ScrollToTop from 'lib/hooks/scroll-to-top';
import { ExtLink } from 'lib/nav';

export const DataPage = () => (
  <article>
    <ScrollToTop />
    <p>
      <Alert severity="info">
        The systemic risk analysis results shown in this tool contain licensed data that must not be shared outside the
        Government of Jamaica. By accessing the tool, you acknowledge that you understand this and agree not to download
        any data or share your access credentials with anyone else.
      </Alert>
    </p>
    <p>
      <Alert
        severity="success"
        action={
          <Button color="inherit" size="small">
            <ExtLink href="https://github.com/nismod/infra-risk-vis/issues">REPORT</ExtLink>
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
            <TableCell>Transport</TableCell>
            <TableCell>Road links and railway lines, ports and airports</TableCell>
            <TableCell>Cost of rehabilitation/reinstating damaged assets</TableCell>
            <TableCell>Rerouting costs + wider effects of service disruption</TableCell>
          </TableRow>
          <TableRow>
            <TableCell>Energy</TableCell>
            <TableCell>Electricity transmission and distribution: generation, lines, poles and substations</TableCell>
            <TableCell>Cost of rehabilitation/reinstating damaged assets</TableCell>
            <TableCell>Wider effects of service disruption</TableCell>
          </TableRow>
          <TableRow>
            <TableCell>Water</TableCell>
            <TableCell>Water supply and wastewater networks, wells and irrigation canals</TableCell>
            <TableCell>Cost of rehabilitation/reinstating damaged assets</TableCell>
            <TableCell>Wider effects of service disruption</TableCell>
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
        <ExtLink href="https://github.com/nismod/infra-risk-vis">github.com/nismod/infra-risk-vis</ExtLink>
      </li>
    </ul>

    <p>The analytics for Jamaica are produced using the code and models here:</p>

    <ul>
      <li>
        <ExtLink href="https://github.com/nismod/jamaica-infrastructure">
          github.com/nismod/jamaica-infrastructure
        </ExtLink>
      </li>
    </ul>
    <h1>Data Sources and Access</h1>

    <p>
      Data comes from multiple sources, including Government of Jamaica bodies, private sector entities, and open data
      sources.
    </p>

    <Alert severity="info">
      The systemic risk analysis results shown in this tool contain licensed data that must not be shared outside the
      Government of Jamaica. By accessing the tool, you acknowledge that you understand this and agree not to download
      any data or share your access credentials with anyone else.
    </Alert>

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
            <TableCell>Fluvial (river) and pluvial (surface) flooding</TableCell>
            <TableCell>
              <img src="/logo-jba.png" alt="JBA Risk Management. The Flood People (R)." width="180" />
              <br />
              <ExtLink href="https://www.jbarisk.com/flood-services/maps-and-analytics/global-flood-maps/">
                JBA global flood map product
              </ExtLink>
            </TableCell>
            <TableCell>1/20, 1/50, 1/100, 1/200, 1/500, and 1/1500</TableCell>
            <TableCell>Flood depths in meters over 30m grid squares</TableCell>
            <TableCell>
              RCP 2.6, 4.5 &amp; 8.5 emission&nbsp;scenarios.
              <br />
              Current + future maps in 2050 and 2080
            </TableCell>
          </TableRow>
          <TableRow>
            <TableCell>Coastal flooding (storm surge)</TableCell>
            <TableCell>
              <ExtLink href="http://web.archive.org/web/20200418050005/https://www.deltares.nl/en/news/regional-study-assess-effects-sea-level-rise-resilience-caribbean/">
                Deltares NL Caribbean product
              </ExtLink>
            </TableCell>
            <TableCell>1/1, 1/2, 1/5, 1/10, 1/50, 1/100</TableCell>
            <TableCell>Flood depths in meters over 90m grid squares </TableCell>
            <TableCell>
              RCP 2.6, 4.5 &amp; 8.5 emission scenarios.
              <br />
              Current + future maps in 2030, 2050, 2070 and 2100
            </TableCell>
          </TableRow>
          <TableRow>
            <TableCell>Tropical cyclones (winds)</TableCell>
            <TableCell>
              <ExtLink href="https://data.4tu.nl/articles/dataset/STORM_IBTrACS_present_climate_synthetic_tropical_cyclone_tracks/12706085/2">
                STORM IBTrACS model
              </ExtLink>
            </TableCell>
            <TableCell>26 different exceedance probabilities from 1/1 to 1/10000</TableCell>
            <TableCell>10 minute sustained maximum wind speeds in m/s at 10km grid squares</TableCell>
            <TableCell>
              RCP 4.5 &amp; 8.5 emission scenarios.
              <br />
              Current + future maps in 2050 and 2100
            </TableCell>
          </TableRow>
          <TableRow>
            <TableCell>Droughts</TableCell>
            <TableCell>REGCM4</TableCell>
            <TableCell>No exceedance probabilities</TableCell>
            <TableCell>Daily rainfall and precipitation</TableCell>
            <TableCell>RCP 2.6, 4.5 &amp; 8.5 emission scenarios</TableCell>
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
            <TableCell>Important failure attributes</TableCell>
            <TableCell>Data sources</TableCell>
          </TableRow>
        </TableHead>
        <TableBody>
          <TableRow>
            <TableCell rowSpan={2}>Energy</TableCell>
            <TableCell>Generation</TableCell>
            <TableCell>9 Power plants</TableCell>
            <TableCell>Damage costs (J$), Population served, GDP disrupted (J$/day)</TableCell>
            <TableCell rowSpan={2}>NSDMD, JPS, MSET, OUR, OpenStreetMap, STATIN</TableCell>
          </TableRow>
          <TableRow>
            <TableCell>Transmission &amp; Distribution</TableCell>
            <TableCell>59 Substations, 30,000 Poles, 11,440 kms of overhead lines</TableCell>
            <TableCell>Damage costs (J$ or J$/m), Population served, GDP disrupted (J$/day)</TableCell>
          </TableRow>
          <TableRow>
            <TableCell rowSpan={4}>Transport</TableCell>
            <TableCell>Airports</TableCell>
            <TableCell>7 airports areas</TableCell>
            <TableCell>Damage costs (J$/m2), Annual passengers, Annual freight (tonnes)</TableCell>
            <TableCell rowSpan={4}>NSDMD, NWA, NROCC, MTM, STATIN</TableCell>
          </TableRow>
          <TableRow>
            <TableCell>Ports</TableCell>
            <TableCell>13 Port dock areas</TableCell>
            <TableCell>Damage costs (J$/m2), Annual passengers, Annual freight (tonnes)</TableCell>
          </TableRow>
          <TableRow>
            <TableCell>Railways</TableCell>
            <TableCell>20 functional stations, 201 kms of functional tracks</TableCell>
            <TableCell>Damage costs (J$ or J$/m), Trade flow disruptions (J$/day)</TableCell>
          </TableRow>
          <TableRow>
            <TableCell>Roads</TableCell>
            <TableCell>572 bridges, 23,200 kms of roads</TableCell>
            <TableCell>
              Damage costs (J$ or J$/m), Reopening costs (J$ or J$/m), Road traffic counts, Trade flow disruptions
              (J$/day)
            </TableCell>
          </TableRow>
          <TableRow>
            <TableCell rowSpan={3}>Water</TableCell>
            <TableCell>Potable water</TableCell>
            <TableCell>1,208 point assets, 10,500 kms of pipelines</TableCell>
            <TableCell>Damage costs (J$ or J$/m), Population served, GDP disrupted (J$/day)</TableCell>
            <TableCell rowSpan={3}>NSDMD, WRA, NWC, NIC, STATIN</TableCell>
          </TableRow>
          <TableRow>
            <TableCell>Irrigation</TableCell>
            <TableCell>178 Wells, 248 kms of canals and 220 kms of pipelines</TableCell>
            <TableCell>Damage costs (J$ or J$/m), Agriculture GDP disrupted (J$/day)</TableCell>
          </TableRow>
          <TableRow>
            <TableCell>Wastewater</TableCell>
            <TableCell>151 Point assets</TableCell>
            <TableCell>Damage costs (J$ or J$/m)</TableCell>
          </TableRow>
          <TableRow>
            <TableCell>Buildings</TableCell>
            <TableCell>
              Commercial, Industrial, Institutional, Mixed Use, Other, Recreation, Residential, Resort
            </TableCell>
            <TableCell>996,682 buildings</TableCell>
            <TableCell>Damage costs (J$/m2), GDP disrupted (J$/day)</TableCell>
            <TableCell>OpenStreetMap, NLA, STATIN</TableCell>
          </TableRow>
        </TableBody>
      </Table>
    </TableContainer>

    <h2>Contextual Map Data</h2>

    <p>
      Background map data is &copy; <ExtLink href="https://www.openstreetmap.org/copyright">OpenStreetMap</ExtLink>{' '}
      contributors, style &copy; <ExtLink href="https://carto.com/attributions">CARTO</ExtLink>.
    </p>

    <p>
      Satellite imagery background is derived from{' '}
      <ExtLink href="https://s2maps.eu">Sentinel-2 cloudless - https://s2maps.eu</ExtLink> by{' '}
      <ExtLink href="https://eox.at">EOX IT Services GmbH</ExtLink> (Contains modified Copernicus Sentinel data 2020).
    </p>
  </article>
);
