import React from 'react';
import PropTypes from 'prop-types';

const Help = ({ topic }) => {
  if (topic === 'flood') {
    return <FloodHelp />;
  }
  if (topic === 'vietnam') {
    return <VietnamHelp />;
  }
  if (topic === 'hazards') {
    return <HazardHelp />;
  }
  return null;
};

Help.propTypes = {
  topic: PropTypes.string.isRequired,
};

const HazardHelp = () => (
  <div className="custom-map-control top-right selected-feature">
    <h4 className="h5">Hazard Data</h4>

    <p>
      The flood (fluvial and coastal) maps show 1/100-year return period flood outlines at 30 arc second grid squares.
    </p>

    <p>The cyclone map shows 1/100-year return period wind speeds at in m/s at 0.1-degree grid squares.</p>

    <p>These are one of several maps used in this study the details of which are described below.</p>

    <table className="table table-sm">
      <thead>
        <tr>
          <th>Hazard</th>
          <th>Probabilities</th>
          <th>Intensities and spatial extents </th>
          <th>Climate scenario information </th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <td>
            Fluvial (river) flooding, from{' '}
            <a target="_blank" rel="noopener noreferrer" href="https://www.wri.org/data/aqueduct-floods-hazard-maps">
              WRI Aqueduct
            </a>
          </td>
          <td rowSpan="2">1/2, 1/5, 1/10, 1/25, 1/50, 1/100, 1/250, 1/500, and 1/1000 </td>
          <td rowSpan="2">Flood depths in meters over 30 arc second grid squares. </td>
          <td rowSpan="2">
            <ul>
              <li>1 current and 5 future climate models</li>
              <li>RCP 4.5 and 8.5 emission scenarios</li>
              <li>Current and future maps in 2030, 2050, 2080</li>
            </ul>
          </td>
        </tr>
        <tr>
          <td>
            Coastal flooding with subsidence (median value), from{' '}
            <a target="_blank" rel="noopener noreferrer" href="https://www.wri.org/data/aqueduct-floods-hazard-maps">
              WRI Aqueduct
            </a>
          </td>
        </tr>
        <tr>
          <td>
            Cyclones, from{' '}
            <a
              target="_blank"
              rel="noopener noreferrer"
              href="https://data.4tu.nl/articles/dataset/STORM_tropical_cyclone_wind_speed_return_periods/12705164/2"
            >
              STORM IBTrACS model
            </a>
          </td>
          <td>28 different probabilities from 1/10 to 1/10000</td>
          <td>3-hour time step wind gust speeds in m/s at 0.1-degree grid squares.</td>
          <td>None</td>
        </tr>
      </tbody>
    </table>
  </div>
);

const FloodHelp = () => (
  <div className="custom-map-control top-right selected-feature">
    <h4 className="h5">Flood Climate Outlooks - Explanation</h4>
    <dl>
      <dt>Baseline</dt>

      <dd>The estimated flooded depths and areas averaged over 1960-1999 based on historical rainfall records.</dd>

      <dt>RCP 4.5</dt>

      <dd>
        The estimated flooded depths and areas averaged over 2010-2049 based on global climate model outputs assuming
        global carbon emissions peak by 2040 before declining. The layer displayed here shows the flood outlines from
        the UK Met Office Hadley Centre Global Environment Model version 2 (HadGEM2-ES) model output, which is 1 of 5
        models used in this study.
      </dd>

      <dt>RCP 8.5</dt>

      <dd>
        The estimated flooded depths and areas averaged over 2010-2049 based on global climate model outputs assuming
        global carbon emissions continue to rise throughout the 21st century. The layer displayed here shows the flood
        outlines from the UK Met Office Hadley Centre Global Environment Model version 2 (HadGEM2-ES) model output,
        which is 1 of 5 models used in this study.
      </dd>
    </dl>
  </div>
);

const VietnamHelp = () => (
  <div className="custom-map-control top-right selected-feature">
    <h4 className="h5">Transport risk and adaptation analysis – Vietnam case-study</h4>

    <p>
      The map here shows the results from a study on transport risk and adaptation analysis of national roads in
      Vietnam, using more detailed information of asset conditions, freight flows, and macroeconomic impacts. The main
      objective of the Vietnam study was to provide a methodological framework to analyse network criticality and
      vulnerability, and to prioritize investments to enhance resilience.
    </p>

    <p>
      Details of the study are here: Oh, J.E., Espinet Alegre, X., Pant, R., Koks, E.E., Russell, T., Schoenmakers, R.
      and Hall, J.W. (2019). Addressing Climate Change in Transport: Volume 2: Pathway to Resilient Transport. World
      Bank, Washington DC. DOI:{' '}
      <a target="_blank" rel="noopener noreferrer" href="http://dx.doi.org/10.1596/32412">
        10.1596/32412
      </a>
      .
    </p>

    <p>
      The map shows the criticality of national roads, highlighted by the line thickness, in terms of the expected
      annual risks (expected annual damages + expected annual economic losses for 30-day disruption durations) due to
      exposures to current fluvial flooding of different return periods and flood depths. Fluvial flood maps were the
      same product as used in the SE Asia study.
    </p>

    <p>The risks were also estimated for future RCP4.5 and RCP8.5 climate change scenarios in 2050.</p>

    <p>By clicking on a road at risk you can also see:</p>

    <p>Flood exposure statistics (across all return periods and for different climate scenarios):</p>

    <ul>
      <li>The range of flood depths at the location of the road</li>
      <li>The minimum return period of exposure for the road, </li>
      <li>The range of road lengths exposed to flooding. </li>
    </ul>

    <p>Criticality metrics (economic impacts of freight flows due to road disruption):</p>

    <ul>
      <li>Increase in rerouting costs or macroeconomic losses in US$/day</li>
      <li>Total economic impact in US$/day</li>
    </ul>

    <p>Risk estimates (for different climate scenarios):</p>

    <ul>
      <li>Expected annual damages to the road link in US$ </li>
      <li>Expected annual economic losses due to freight disruptions along the link in US$/day</li>
    </ul>

    <p>Adaptation option:</p>

    <ul>
      <li>The option to upgrade to a climate resilience road that eliminates all risks</li>
    </ul>

    <p>Adaptation cost estimates (across all return periods and for different climate scenarios over 35 years):</p>

    <ul>
      <li>Ranges of initial investment in US$ required for the road a href make it climate resilient</li>
      <li>Ranges of net present value of maintenance costs in US$ require for the road link over the years</li>
      <li>Ranges of net present value total adaptation investment in US$ for the road link over the years </li>
    </ul>

    <p>Adaptation cost estimates per km (for different climate scenarios over 35 years):</p>

    <ul>
      <li>The adaptation cost estimate divided by the length of the road links</li>
    </ul>

    <p>
      Benefit-cost ratio estimates (for different climate scenarios over 35 years, and sensitivity to durations of
      disruptions and future GDP growth scenarios):
    </p>

    <ul>
      <li>
        Cost – Net present value of maximum total adaptation investment in US$ for the road link over the years needed
        for eliminating climate risks
      </li>
      <li>Benefit – Net present value of range of avoided risks in US$ for the road link over the years</li>
      <li>Benefit-Cost Ratio – Benefit divided by the Cost</li>
    </ul>
  </div>
);

export default Help;
