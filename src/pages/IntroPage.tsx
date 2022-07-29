import React from 'react';
import { styled } from '@mui/material/styles';
import { Divider, Grid, Paper, Stack, Typography } from '@mui/material';
import ScrollToTop from 'lib/hooks/scroll-to-top';

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
  <article>
    <ScrollToTop />
    <Grid container columnSpacing={8} rowSpacing={4}>
      <Grid item xs={6}>
        <HeadingBox sx={{mt:-2, pt:8}}>
          <Typography variant="h1">

            Climate-related risk analytics for transport, energy &amp; water
            infrastructure in the Caribbean

          </Typography>
        </HeadingBox>
      </Grid>
      <Grid item xs={6}>
        <TextBox sx={{mt:-2, pt:8}}>
          <p>
            The Systemic Risk Assessment Tool (SRAT) supports
            climate adaptation decision-making by identifying spatial
            criticalities and risks under current and future climate scenarios.
          </p>
          <Typography variant="h2">Transport</Typography>
          <p>Roads, from major to tertiary.</p>
          <Typography variant="h2">Energy</Typography>
          <p>Electricity generation and transmission.</p>

        </TextBox>
      </Grid>
      <Grid item xs={12}>
        <TextBox sx={{backgroundColor: 'rgba(255, 255, 255, 0.9)'}}>

          <p>The research, analysis and development has been led by
            researchers in the&nbsp;<a href="https://opsis.eci.ox.ac.uk/" target="_blank" rel="noopener noreferrer">Oxford Programme for
            Sustainable Infrastructure Systems</a>, University of Oxford,
            funded by UK Aid through the UK Foreign and
            Commonwealth Development Office (FCDO)</p>
          <p>
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
            <a href="https://www.gov.uk/guidance/uk-aid" target="_blank" rel="noopener noreferrer">
              <img height="100" src="/logo-ukaid.png" alt="UK AID" />
            </a>
          </Stack>
          </p>

          <p><small>
            Photo credit: Hurricane Irma, 7 September 2017.
            Data: MODIS/Terra (NASA WorldView). Processed by Antti Lipponen&nbsp;(
            <a href="https://twitter.com/anttilip" target="_blank" rel="noopener noreferrer">@anttilip</a>)&nbsp;
            <a href="https://creativecommons.org/licenses/by/2.0/" target="_blank" rel="noopener noreferrer">CC-BY</a>
          </small></p>

        </TextBox>

      </Grid>
    </Grid>
  </article>
);
