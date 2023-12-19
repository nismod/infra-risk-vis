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

      <Typography variant="h2">Technical Reports</Typography>

      <p>
        This is the decision support tool produced as part of project HVT043, “Decision Support Systems for Resilient
        Strategic Transport Networks in Low Income Countries”. It is built around an interactive web platform and aims
        to support investment decisions and options selection for long distance strategic land transport networks
        exposed to climate risks. It is the first multi-state transport infrastructure decision support system in a
        low-income country context, based on a case study region covering Uganda, Zambia, Kenya and Tanzania.
      </p>
      <p>
        The summary and final reports provide an overview of the research findings underpinning the decision support
        tool which has been developed during the project. The underlying research has focused on developing a range of
        future background scenarios for transport development in the case study region, identifying and assembling
        datasets which form the basis for an assessment of transport resilience and sustainability.
      </p>
      <p>
        Data requirements, methodologies, related frameworks and example results for the underlying research are
        presented throughout the report, which also summarises the development of the decision support tool and provides
        a case study example based on a road enhancement project in Kenya. The case study was one of five identified in
        discussions during stakeholder workshops. Details of the three sets of online workshops are provided, as well as
        an overview of the four in-country demonstration workshops carried out in September 2022.
      </p>
      <ul>
        <li>
          Summary report{' '}
          <a
            href="https://transport-links.com/hvt-publications/summary-report-decision-support-systems-for-resilient-strategic-transport-networks-in-low-income-countries"
            target="_blank"
            rel="noopener noreferrer"
          >
            [external]
          </a>
          ,{' '}
          <a href="/summary-report-decision-support-systems-for-resilient-strategic-transport-networks-in-low-income-countries.pdf">
            [pdf]
          </a>
        </li>

        <li>
          Final report{' '}
          <a
            href="https://transport-links.com/hvt-publications/final-report-decision-support-systems-for-resilient-strategic-transport-networks-in-low-income-countries"
            target="_blank"
            rel="noopener noreferrer"
          >
            [external]
          </a>
          ,{' '}
          <a href="/final-report-decision-support-systems-for-resilient-strategic-transport-networks-in-low-income-countries.pdf">
            [pdf]
          </a>
        </li>
        <li>
          User guide{' '}
          <a
            href="https://transport-links.com/hvt-publications/sustainability-and-risk-analysis-tool-user-guide"
            target="_blank"
            rel="noopener noreferrer"
          >
            [external]
          </a>
          , <a href="sustainability-and-risk-analysis-tool-user-guide.pdf">[pdf]</a>
        </li>
      </ul>

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
