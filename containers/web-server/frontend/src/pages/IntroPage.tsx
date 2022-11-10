import { Alert, Divider, Grid, Paper, Stack, Typography, useMediaQuery } from '@mui/material';

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

export const IntroPage = () => {
  const isMobile = useMediaQuery((theme: any) => theme.breakpoints.down('md'));

  return (
  <div className="home">
    <article>
      <ScrollToTop />
      <Grid container columnSpacing={8} rowSpacing={4}>
        <Grid item md={6} sx={{width:'100%'}}>
          <HeadingBox sx={{ mt: -2, pt: 8 }}>
            <Typography variant="h1">Global climate-related risk analytics</Typography>
          </HeadingBox>
        </Grid>
        <Grid item md={6}>
          <TextBox sx={{ mt: -2, py: 8 }}>
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
        {
          isMobile?
          <Grid item xs={12}>
          <Alert sx={{my:2}} severity='warning'>

            This site is not currently well-designed for small screens. For
            a better experience, we recommend visiting from a larger device
            or window.

          </Alert>
          </Grid>
          : null
        }
        <Grid item xs={12}>
          <TextBox sx={{ backgroundColor: 'rgba(255, 255, 255, 0.9)' }}>
            <p>

              The research, analysis and development has been led by researchers
              in the <a href="https://opsis.eci.ox.ac.uk/" target="_blank"
                rel="noopener noreferrer"> Oxford Programme for Sustainable
                Infrastructure Systems</a> at the University of Oxford as part of
              the <a href="http://www.globalresilienceindex.org/"
                target="_blank" rel="noopener noreferrer">Global Resilience Index
                Initiative</a>.

            </p>

            <Stack
              direction={{ xs: 'column', md: 'row' }}
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


              <a href="http://www.globalresilienceindex.org/" target="_blank" rel="noopener noreferrer">
                <img
                  height="80"
                  src="/logo-grii.png"
                  alt="Global Resilience Index Initiative"
                />
              </a>
            </Stack>

            <Typography variant="h2">GRII Partners</Typography>


            <Stack
              direction={{ xs: 'column', md: 'row' }}
              divider={<Divider orientation="vertical" flexItem />}
              spacing={{ md: 1, lg: 2 }}
              sx={{ mb: 2 }}
              justifyContent="center"
              alignItems="center"
            >
              <a href="https://ukcgfi.org/" target="_blank" rel="noopener noreferrer">
                <img height="50" src="/logo-cgfi.png" alt="CGFI" />
              </a>
              <a href="https://resilientinvestment.org/" target="_blank" rel="noopener noreferrer">
                <img height="70" src="/logo-ccri.png" alt="CCRI" />
              </a>
              <a href="https://www.cdri.world/" target="_blank" rel="noopener noreferrer">
                <img height="60" src="/logo-cdri.png" alt="CDRI" />
              </a>
              <a href="https://www.globalquakemodel.org/who-we-are" target="_blank" rel="noopener noreferrer">
                <img height="50" src="/logo-gem.png" alt="GEM" />
              </a>
              <a href="https://www.insdevforum.org/" target="_blank" rel="noopener noreferrer">
                <img height="30" src="/logo-idf.png" alt="IDF" />
              </a>
              <a href="https://www.undrr.org/" target="_blank" rel="noopener noreferrer">
                <img height="40" src="/logo-undrr.png" alt="UNDRR" />
              </a>
            </Stack>

            <Typography variant="h2">Funding and support</Typography>

            <p>

            This project is led by researchers in the <a
            href="https://opsis.eci.ox.ac.uk/" target="_blank" rel="noopener
            noreferrer">Oxford Programme for Sustainable Infrastructure
            Systems</a> in the Environmental Change Institute, University of
            Oxford, with contributions of data from the <a
            href="https://www.globalquakemodel.org/who-we-are" target="_blank"
            rel="noopener noreferrer">Global Earthquake Model Foundation</a>,
            the <a href="https://www.cgfi.ac.uk/spatial-finance-initiative/"
            target="_blank" rel="noopener noreferrer">Spatial Finance
            Initiative</a>, as well as the many open data sources listed <Link
            to="/data">here</Link>. Appropriate open data sources are regularly
            reviewed as part of the GRII taskforce.

            </p>
            <p>

            This project has been funded by the UK Natural Environment Research
            Council (NERC) through the UK Centre for Greening Finance and
            Investment, the World Bank Group, Insurance for Development Forum,
            and Willis Towers Watson.

            </p>

            <Stack
              direction={{ xs: 'column', md: 'row' }}
              divider={<Divider orientation="vertical" flexItem />}
              spacing={2}
              sx={{ mb: 2 }}
              justifyContent="center"
              alignItems="center"
            >
              <a href="https://www.ukri.org/councils/nerc/" target="_blank" rel="noopener noreferrer">
                <img height="60" src="/logo-nerc.png" alt="NERC" />
              </a>
              <a href="https://www.worldbank.org/en/home" target="_blank" rel="noopener noreferrer">
                <img height="40" src="/logo-wbg.png" alt="World Bank Group" />
              </a>
              <a href="https://www.insdevforum.org/" target="_blank" rel="noopener noreferrer">
                <img height="40" src="/logo-idf.png" alt="IDF" />
              </a>
              <a href="https://www.wtwco.com/" target="_blank" rel="noopener noreferrer">
                <img height="30" src="/logo-wtw.png" alt="WTW" />
              </a>
            </Stack>

            <p>

            This builds on previous research and development funded by UK AID
            through the UK Foreign and Commonwealth Development Office (FCDO) as
            part of a project with the Government of Jamaica (GoJ) under the
            Coalition for Climate Resilient Investment&rsquo;s (CCRI) work on
            "Systemic Resilience" in collaboration with the Green Climate Fund,
            and also through the High-Volume Transport Applied Research project.

            </p>
            <p>

            Similarly, earlier versions of the tool piloted in Argentina and
            South-East Asia were funded by the Disaster Risk Financing and
            Insurance Program (DRFIP) of the World Bank with support from the
            Japan—World Bank Program for Mainstreaming DRM in Developing
            Countries, which is financed by the Government of Japan and managed
            by the Global Facility for Disaster Reduction and Recovery (GFDRR)
            through the Tokyo Disaster Risk Management Hub.

            </p>

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
  )
};
