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
  <div className="home">
    <article>
      <ScrollToTop />
      <Grid container columnSpacing={8} rowSpacing={4}>
        <Grid item xs={6}>
          <HeadingBox sx={{ mt: -2, pt: 8 }}>
            <Typography variant="h1">

              Global Climate-related risk analytics
            </Typography>
          </HeadingBox>
        </Grid>
        <Grid item xs={6}>
          <TextBox sx={{ mt: -2, pt: 8 }}>
            <p>
              The Global Systemic Risk Assessment Tool (G&#8209;SRAT) supports
              climate adaptation decision-making by identifying spatial
              criticalities and risks under current and future climate scenarios.
            </p>
            <Typography variant="h2">Transport</Typography>
            <p>Roads, rail, ports and airports.</p>
            <Typography variant="h2">Energy</Typography>
            <p>Electricity generation, transmission and distribution.</p>
            <Typography variant="h2">Water</Typography>
            <p>Water supply, wastewater and irrigation.</p>

          </TextBox>
        </Grid>
        <Grid item xs={12}>
          <TextBox sx={{ backgroundColor: 'rgba(255, 255, 255, 0.9)' }}>

            <p>The research, analysis and development has been led by
              researchers in the&nbsp;<a href="https://opsis.eci.ox.ac.uk/" target="_blank" rel="noopener noreferrer">Oxford Programme for
                Sustainable Infrastructure Systems</a>, University of Oxford</p>
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
              </Stack>
            </p>


            <Typography variant="h2">Funding and support</Typography>

            <p>This project is led by researchers in the&nbsp;
              <a href="https://opsis.eci.ox.ac.uk/" target="_blank" rel="noopener noreferrer">Oxford Programme for
                Sustainable Infrastructure Systems</a> in the Environmental Change
              Institute, University of Oxford.</p>

            <Stack
              direction="row"
              divider={<Divider orientation="vertical" flexItem />}
              spacing={2}
              justifyContent="center"
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

          </TextBox>

        </Grid>
      </Grid>
    </article>
  </div>
);
