import React from 'react';
import { Alert, Button, Divider, Paper, Stack, Table, TableBody, TableCell, TableContainer, TableHead, TableRow } from '@mui/material';
import { Link } from 'react-router-dom';
import ScrollToTop from 'lib/hooks/scroll-to-top';

export const IntroPage = () => (
  <article>
    <ScrollToTop />
    <h1>Caribbean Systemic Risk Assessment Tool</h1>

    <p>The Systemic Risk Assessment Tool (SRAT) presents
    climate-related risk analytics for transport and energy
    infrastructure across the Caribbean.</p>

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

    <p>For more details on the infrastructure and hazard data used, see
    the <Link to="/data">Data page</Link></p>

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

    <h2>Funding and support</h2>

    <p>This research is led by researchers in the Oxford Programme for
    Sustainable Infrastructure Systems in the Environmental Change Institute,
    University of Oxford as part of a project funded by UK Aid (FCDO).</p>

    <Stack
      direction="row"
      divider={<Divider orientation="vertical" flexItem />}
      spacing={2}
      justifyContent="center"
      alignItems="center"
    >
      <a href="https://opsis.eci.ox.ac.uk" target="_blank" rel="noopener noreferrer">
        <img height="100" src="/logo-opsis.png" alt="Oxford Programme for Sustainable Infrastructure Systems" />
      </a>
      <a href="https://www.gov.uk/guidance/uk-aid" target="_blank" rel="noopener noreferrer">
        <img height="100" src="/logo-ukaid.png" alt="UK AID" />
      </a>
    </Stack>
  </article>
);
