import React from 'react';
import { Alert, Button, Divider, Paper, Stack, Table, TableBody, TableCell, TableContainer, TableHead, TableRow } from '@mui/material';
import { Link } from 'react-router-dom';
import ScrollToTop from 'lib/hooks/scroll-to-top';

export const IntroPage = () => (
  <article>
    <ScrollToTop />
    <h1>Global Systemic Risk Assessment Tool</h1>

    <p>The Global Systemic Risk Assessment Tool (G-SRAT) presents
      climate-related risk analytics for transport, energy and water
      infrastructure globally.</p>

    <Alert severity="info">TODO: Licensing</Alert>

    <br></br>

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

    <p>For more details on the infrastructure and hazard data used, see
      the <Link to="/data">Data page</Link></p>

    <p>The primary output metrics from the analysis are:</p>

    <ul>
      <li>Expected Annual Damages (EAD) (direct physical risks) estimated as the
        area under the direct damage vs exceedance probability curve </li>
      <li>Expected Annual Economic Losses (EAEL) (indirect economic risks)
        estimated as the area under the economic loss vs exceedance probability
        curve </li>
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

    <p>The global analytics are produced using the code and models here:</p>

    <ul>
      <li>
        <a href="TODO" target="blank">
          TODO
        </a>
      </li>
    </ul>

    <h2>Collaboration and Engagement</h2>

    <p>TODO</p>

    <Stack
      direction="row"
      divider={<Divider orientation="vertical" flexItem />}

      justifyContent="center"
      alignItems="center"
      spacing={2}
    >
      <a href="https://opsis.eci.ox.ac.uk" target="_blank" rel="noopener noreferrer">
        <img height="100" src="/logo-opsis.png" alt="Oxford Programme for Sustainable Infrastructure Systems" />
      </a>
    </Stack>


    <h2>Funding and support</h2>

    <p>TODO</p>

    <p>The initiative forms part of the Coalition for Climate Resilient
      Investment&rsquo;s (CCRI) and collaboration with the Green Climate Fund.</p>

    <Stack
      direction="row"
      divider={<Divider orientation="vertical" flexItem />}
      spacing={2}
      justifyContent="right"
      alignItems="center"
    >
      <a href="https://www.gov.uk/guidance/uk-aid" target="_blank" rel="noopener noreferrer">
        <img height="100" src="/logo-ukaid.png" alt="UK AID" />
      </a>
      <a href="https://www.greenclimate.fund/" target="_blank" rel="noopener noreferrer">
        <img height="100" src="/logo-gcf.png" alt="Green Climate Fund" />
      </a>
      <a href="https://resilientinvestment.org/" target="_blank" rel="noopener noreferrer">
        <img height="100" src="/logo-ccri.png" alt="Coalition for Climate Resilient Investment" />
      </a>
    </Stack>
  </article>
);
