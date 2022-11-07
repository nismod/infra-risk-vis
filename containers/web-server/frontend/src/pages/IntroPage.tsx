import { Divider, Grid, Paper, Stack, Typography } from '@mui/material';
import { styled } from '@mui/material/styles';
import React from 'react';

import ScrollToTop from '@/lib/react/ScrollToTop';
import { Link } from 'react-router-dom';

const HeadingBox = styled(Paper)(({ theme }) => ({
  backgroundColor: 'rgba(0, 92, 97, 0.3)',
  color: '#fff',
  padding: theme.spacing(2),
  borderRadius: 0,
}));

const TextBox = styled(Paper)(() => ({
  backgroundColor: 'rgba(194, 219, 231, 0.9)',
  color: '#333',
  padding: '16px 32px',
  borderRadius: 0,
}));

export const IntroPage = () => (
  <div className="home">
    <article>
      <ScrollToTop />
      <Grid container columnSpacing={8} rowSpacing={4}>
        <Grid item xs={6}>
          <HeadingBox sx={{ mt: -2, pt: 8 }}>
            <Typography variant="h1">Global climate-related risk analytics</Typography>
          </HeadingBox>
        </Grid>
        <Grid item xs={6}>
          <TextBox sx={{mt:-2, pt:8}}>
            <p>

              The Global Systemic Risk Assessment Tool (G-SRAT)
              is the Data and Analytics Portal for the <a
              href="http://www.globalresilienceindex.org/">Global Resilience Index Initative
              (GRII)</a>.

            </p>
            <p>
              This tool aims to support
              climate adaptation decision-making by identifying spatial
              vulnerabilities and risks under current and future climate scenarios.
            </p>

          </TextBox>
        </Grid>
        <Grid item xs={12}>
          <TextBox sx={{ backgroundColor: 'rgba(255, 255, 255, 0.9)' }}>
            <p>
              The research, analysis and development has been led by researchers in the&nbsp;
              <a href="https://opsis.eci.ox.ac.uk/" target="_blank" rel="noopener noreferrer">
                Oxford Programme for Sustainable Infrastructure Systems
              </a>
              , University of Oxford
            </p>
            <p>
              <Stack
                direction="row"
                divider={<Divider orientation="vertical" flexItem />}
                justifyContent="center"
                alignItems="center"
                spacing={2}
              >
                <a href="https://opsis.eci.ox.ac.uk" target="_blank" rel="noopener noreferrer">
                  <img
                    height="100"
                    src="/logo-opsis.png"
                    alt="Oxford Programme for Sustainable Infrastructure Systems"
                  />
                </a>
              </Stack>
            </p>

            <Typography variant="h2">Funding and support</Typography>

            <p>
              This project is led by researchers in the&nbsp;
              <a href="https://opsis.eci.ox.ac.uk/" target="_blank" rel="noopener noreferrer">
                Oxford Programme for Sustainable Infrastructure Systems
              </a>{' '}
              in the Environmental Change Institute, University of Oxford.
            </p>

            <Stack
              direction="row"
              divider={<Divider orientation="vertical" flexItem />}
              spacing={2}
              justifyContent="center"
              alignItems="center"
            ></Stack>

            <Typography variant="h2">Disclaimer</Typography>

            <p>This tool is provided for general information only and is not
            intended to amount to advice on which you should rely. You must
            obtain professional or specialist advice before taking, or
            refraining from, any action on the basis of the content on our site.
            </p>

            <p> Although we make reasonable efforts to update the information on
            our site, we make no representations, warranties or guarantees,
            whether express or implied, that the content on our site (including
            this tool) is accurate, complete or up to date. The University of
            Oxford accepts no liability in relation to any issues or liabilities
            that may subsequently arise from use of the data or this tool for
            any purpose. Please consult our <Link to="/terms-of-use">website
            terms of use</Link> for more information about our liability to you.
            </p>

          </TextBox>
        </Grid>
      </Grid>
    </article>
  </div>
);
