import React from 'react';
import { Alert, Button, Paper, Table, TableBody, TableCell, TableContainer, TableHead, TableRow } from '@mui/material';

export const IntroPage = () => (
  <article>
    <h1>Jamaica Systemic Risk Assessment Tool</h1>

    <p>This prototype tool presents infrastructure risk analytics for
    Jamaica.</p>

    <Alert severity="info">The systemic risk analysis results shown in this
    tool contain licensed data that must not be shared outside the Government of
    Jamaica. By accessing the tool, you acknowledge that you understand this and
    agree not to download any data or share your access credentials with anyone
    else.</Alert>

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

    The tool is under active development. Please tell us if anything is not
    working as it should and suggest potential improvements.

    </Alert>

    <p>The modelling and analysis presented here aim to support climate
    adaptation decision-making by identifying spatial criticalities and risks
    under current and future hazard scenarios. It comprises a direct damage
    estimation and an indirect economic loss estimation of GDP disruptions due
    to asset failures and service disruption.</p>

    <TableContainer component={Paper}>
      <Table aria-label="simple table">
        <TableHead>
          <TableRow>
            <TableCell>Infrastructure</TableCell>
            <TableCell>Assets</TableCell>
            <TableCell>Expected Annual Damages (EAD)</TableCell>
            <TableCell>Expected Annual Economic Losses (EAEL)</TableCell>
          </TableRow>
        </TableHead>
        <TableBody>
          <TableRow>
            <TableCell>Transport</TableCell>
            <TableCell>Road links and railway tracks</TableCell>
            <TableCell>Cost of rehabilitation/reinstating damaged assets</TableCell>
            <TableCell>Rerouting costs + wider effects of service disruption</TableCell>
          </TableRow>
          <TableRow>
            <TableCell>Energy</TableCell>
            <TableCell>Electricity transmission and distribution grid: generation, lines and substations</TableCell>
            <TableCell>Cost of rehabilitation/reinstating damaged assets</TableCell>
            <TableCell>Wider effects of service disruption</TableCell>
          </TableRow>
          <TableRow>
            <TableCell>Water</TableCell>
            <TableCell>Water supply and wastewater networks, abstraction points, irrigation schemes</TableCell>
            <TableCell>Cost of rehabilitation/reinstating damaged assets</TableCell>
            <TableCell>Wider effects of service disruption</TableCell>
          </TableRow>
        </TableBody>
      </Table>
    </TableContainer>

    <p>This tool to visualize the model outputs is developed and documented here:</p>

    <ul>
      <li>
        <a href="https://github.com/nismod/infra-risk-vis" target="blank">
          github.com/nismod/infra-risk-vis
        </a>
      </li>
    </ul>

    <p>The analytics for Jamaica are produced using the code and models here:</p>

    <ul>
      <li>
        <a href="https://github.com/nismod/jamaica-infrastructure" target="blank">
          github.com/nismod/jamaica-infrastructure
        </a>
      </li>
    </ul>

    <h2>Funding support</h2>

    <p>
      This research is led by researchers in the Oxford Programme for Sustainable Infrastructure Systems in the
      Environmental Change Institute, University of Oxford, for the Government of Jamaica (GoJ) as part of a project
      funded by UK Aid (FCDO). The initiative forms part of the Coalition for Climate Resilient Investment’s (CCRI)
      collaboration with the GoJ, which includes analysis of nature-based approaches to build resilience in Jamaica to
      be procured and funded by the Green Climate Fund (GCF).
    </p>
  </article>
);
