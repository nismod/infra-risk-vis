import { Stack, Typography } from '@mui/material';
import ScrollToTop from 'lib/hooks/scroll-to-top';

import { HeadingBox } from './HeadingBox';

export const IntroPage = () => (
  <article>
    <ScrollToTop />
    <HeadingBox>
      <Typography variant="h1">
        Climate-related risk analytics for road &amp; rail transport infrastructure in East Africa
      </Typography>
    </HeadingBox>
    <div className="home" style={{ height: '16rem' }}></div>
    <div className="centred">
      <hr className="minibar" />
      <Typography variant="h5" component="p">
        The Systemic Risk Assessment Tool (SRAT) supports climate adaptation decision-making by identifying spatial
        criticalities and risks under current and future climate scenarios.
      </Typography>
      <Typography variant="h5" component="p">
        This tool presents analysis of major and minor roads, rail lines and stations in four countries: Kenya,
        Tanzania, Uganda and Zambia.
      </Typography>

      <Stack
        direction={{ xs: 'column', sm: 'row' }}
        flexWrap="wrap"
        spacing={2}
        sx={{ mt: 4, mb: 2 }}
        justifyContent="center"
        alignItems="center"
      >
        <a
          href="https://www.southampton.ac.uk/research/groups/transportation-group"
          target="_blank"
          rel="noopener noreferrer"
        >
          <img height="80" src="/logo-southampton.png" alt="University of Southampton Transportation Research Group" />
        </a>
        <a href="https://opsis.eci.ox.ac.uk" target="_blank" rel="noopener noreferrer">
          <img height="110" src="/logo-opsis.png" alt="Oxford Programme for Sustainable Infrastructure Systems" />
        </a>
      </Stack>

      <p>
        The research, analysis and development has been led by researchers in the University of Southampton's{' '}
        <a
          href="https://www.southampton.ac.uk/research/groups/transportation-group"
          target="_blank"
          rel="noopener noreferrer"
        >
          Transportation Research Group
        </a>{' '}
        and the{' '}
        <a href="https://opsis.eci.ox.ac.uk/" target="_blank" rel="noopener noreferrer">
          Oxford Programme for Sustainable Infrastructure Systems
        </a>
        , University of Oxford, supported by engagement with infrastructure and climate specialists and related
        government bodies.
      </p>

      <Typography variant="h2">Funding and support</Typography>

      <p>
        This research was funded by UKAID through the UK Foreign, Commonwealth &amp; Development Office under the High
        Volume Transport Applied Research Programme, which is managed by DT Global.{' '}
      </p>

      <Stack
        direction={{ xs: 'column', sm: 'row' }}
        flexWrap="wrap"
        spacing={2}
        sx={{ mt: 4, mb: 2 }}
        justifyContent="center"
        alignItems="center"
      >
        <a href="https://www.gov.uk/guidance/uk-aid" target="_blank" rel="noopener noreferrer">
          <img height="100" src="/logo-ukaid.png" alt="UK AID" />
        </a>
        <a href="https://transport-links.com/" target="_blank" rel="noopener noreferrer">
          <img height="100" src="/logo-hvt.png" alt="High Volume Transport Applied Research" />
        </a>
        <a href="https://transport-links.com/" target="_blank" rel="noopener noreferrer">
          <img height="60" src="/logo-dt.svg" alt="DT Global" />
        </a>
      </Stack>

      <p>Photo credit: Flooding in the Tana River, Kenya, 1988. US National Archives, Public domain.</p>
    </div>
  </article>
);
