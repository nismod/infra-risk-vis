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
        <HeadingBox sx={{mt:-2, pt:8}}>
          <Typography variant="h1">

            Climate-related risk analytics for transport, energy &amp; water
            infrastructure in East Africa

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
          <Typography variant="h2">Road</Typography>
          <p>Major and minor roads.</p>
          <Typography variant="h2">Rail</Typography>
          <p>Rail lines and stations.</p>

        </TextBox>
      </Grid>
      <Grid item xs={12}>
        <TextBox sx={{backgroundColor: 'rgba(255, 255, 255, 0.9)'}}>

          <p>The research, analysis and development has been led by researchers
          in the University of Southampton's <a
          href="https://www.southampton.ac.uk/research/groups/transportation-group"
          target="_blank" rel="noopener noreferrer">Transportation Research
          Group</a> and the <a href="https://opsis.eci.ox.ac.uk/"
          target="_blank" rel="noopener noreferrer">Oxford Programme for
          Sustainable Infrastructure Systems</a>, University of Oxford,
          supported by engagement with infrastructure and climate specialists
          and related government bodies.</p>

          <Stack
            direction="row"
            divider={<Divider orientation="vertical" flexItem />}

            justifyContent="center"
            alignItems="center"
            spacing={2}
            my={2}
            >
            <a href="https://www.southampton.ac.uk/research/groups/transportation-group"
              target="_blank" rel="noopener noreferrer">
              <img height="80" src="/logo-southampton.png" alt="University of Southampton Transportation Research Group" />
            </a>
            <a href="https://opsis.eci.ox.ac.uk" target="_blank" rel="noopener noreferrer">
              <img height="110" src="/logo-opsis.png" alt="Oxford Programme for Sustainable Infrastructure Systems" />
            </a>
          </Stack>


          <Typography variant="h2">Funding and support</Typography>

          <p>This research was funded by UKAID through the UK Foreign,
          Commonwealth &amp; Development Office under the High Volume Transport
          Applied Research Programme, managed by IMC Worldwide. </p>

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
            <a href="https://www.imcworldwide.com/project/hvt/" target="_blank" rel="noopener noreferrer">
              <img height="100" src="/logo-hvt.png" alt="High Volume Transport Applied Research"/>
            </a>
            <a href="https://www.imcworldwide.com/project/hvt/" target="_blank" rel="noopener noreferrer">
              <img height="90" src="/logo-imc.png" alt="IMC Worldwide"/>
            </a>
          </Stack>

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
  </div>
);
