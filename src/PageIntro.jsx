import React from 'react';

const PageIntro = () => (
  <article>
  <h1>Southeast Asia Infrastructure Risk Prototype</h1>

  <p>This prototype tool presents infrastructure risk analytics for South-East
  Asia (SEA) using open-data sources on fluvial and coastal flooding hazard maps
  along with cyclone hazard maps. The risks are analysed and visualised for
  power plants, electricity transmission lines, road networks, railway networks,
  ports and airports in SEA. The analysis has been undertaken for the World Bank
  by Oxford Infrastructure Analytics Ltd.</p>

  <p>The purpose of the prototype is to illustrate the type and nature of
  simulation results and network data currently available in the SEA region, and
  how &ndash; using the existing toolchain &ndash; this data might be accessed
  and interrogated at a sub-national, national or super-national scale.</p>

  <p>The modelling and analysis presented here aim to support Disaster Risk
  Finance and Insurance (DRFI) decision-making by identifying spatial
  criticalities and risks under current and future hazard scenarios. It
  comprises a direct damage estimation and an indirect economic loss estimation
  of GDP disruptions due to asset failures and service disruption.</p>

  <table className="table table-sm table-striped">
    <thead>
      <th>Infrastructure</th>
      <th>Assets</th>
      <th>Expected Annual Damages (EAD)</th>
      <th>Expected Annual Economic Losses (EAEL)</th>
    </thead>
    <tbody>
      <tr>
        <td>Road</td>
        <td>Road links</td>
        <td>Cost of rehabilitation/reinstating damaged assets</td>
        <td>Rerouting costs* + Macroeconomic losses**</td>
      </tr>
      <tr>
        <td>Rail</td>
        <td>Railway tracks</td>
        <td>Cost of rehabilitation/reinstating damaged assets</td>
        <td>Rerouting costs* + Macroeconomic losses**</td>
      </tr>
      <tr>
        <td>Electricity</td>
        <td>Electricity lines</td>
        <td>Cost of rehabilitation/reinstating damaged assets</td>
        <td>Macroeconomic losses**</td>
      </tr>
    </tbody>
  </table>

  <p><small><em>* Rerouting costs are computed for the Vietnam case study only and
  ignored for other regions. The pages showing results for all of Southeast Asia
  do not include the rerouting costs for Vietnam.</em></small></p>

  <p><small><em>** Macroeconomic losses are computed based on GDP for all regions in the
  Southeast Asia analysis. In the Vietnam case study, IO models were used to
  estimate wider effects of disruption on the macroeconomy - these results are
  shown only on the Vietnam case study page.</em></small></p>

  <p>The concepts and model results presented here are documented in the study
  report:</p>

  <ul><li>Pant, R., Russell, T., Glasgow, G., Verschuur, J., Gavin, H., Fowler,
  T. &amp; Hall, J.W. (2021). <em>Analytics for Financial Risk Management of
  Critical Infrastructure in Southeast Asia – Final Report.</em> Oxford
  Infrastructure Analytics Ltd., Oxford, UK. (Available on request from the
  World Bank)</li></ul>

  <p>The tool being used to visualize the model outputs is developed and
  documented here:</p>

  <ul><li><a href="https://github.com/oi-analytics/oi-risk-vis"
  target="blank">github.com/oi-analytics/oi-risk-vis</a></li></ul>

  <p>The Southeast Asia analytics are produced using the code here:</p>

  <ul><li><a href="https://github.com/oi-analytics/seasia"
  target="blank">github.com/oi-analytics/seasia</a></li></ul>

  <p>The outputs specific to Vietnam visualized here were generated from a model
  created and documented here:</p>

  <ul><li><a href="https://github.com/oi-analytics/vietnam-transport"
  target="blank">github.com/oi-analytics/vietnam-transport</a></li></ul>

  <h2>Funding support</h2>

  <p>This project is led by the Disaster Risk Financing and Insurance Program
  (DRFIP) of the World Bank with support from the Japan&mdash;World Bank Program
  for Mainstreaming DRM in Developing Countries, which is financed by the
  Government of Japan and managed by the Global Facility for Disaster Reduction
  and Recovery (GFDRR) through the Tokyo Disaster Risk Management Hub.</p>

  </article>
);

export default PageIntro;
