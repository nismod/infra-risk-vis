import { Divider, Grid, Paper, Stack, Typography } from '@mui/material';
import { styled } from '@mui/material/styles';
import React from 'react';

import ScrollToTop from '@/lib/react/ScrollToTop';

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
            <Typography variant="h1">Global Climate-related risk analytics</Typography>
          </HeadingBox>
        </Grid>
        <Grid item xs={6}></Grid>
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
          </TextBox>
        </Grid>
      </Grid>
    </article>
  </div>
);
