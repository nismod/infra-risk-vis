import React from 'react';

const PageIntro = () => (
  <article>
  <h1>Southeast Asia Infrastructure Risk: Prototype</h1>

  <p>This prototype tool presents infrastructure risk analytics for South-East
  Asia (SEA) using open-data sources on fluvial and coastal flooding hazard maps
  along with cyclone hazard maps. The risks are analysed and visualised for
  power plants, electricity transmission lines, road networks, railway networks,
  ports and airports in SEA. The analysis has been undertaken for the World Bank
  by Oxford Infrastructure Analytics Ltd.</p>

  <p>The purpose of the prototype is to illustrate the type and nature of
  simulation results and network data currently available in the SEA region, and
  how &ndash; using the existing toolchain &ndash; this data might be accessed
  and interrogated at a sub-national, national or super-national scale. </p>

  <p>The modelling and analysis presented here aim to support Disaster Risk
  Finance and Insurance (DRFI) decision-making by identifying spatial
  criticalities and risks under current and future hazard scenarios. It
  comprises a direct damage estimation and an indirect economic loss estimation
  of GDP disruptions due to asset failures and service disruption.</p>

  <p>The concepts and model results presented here are documented in the study
  report:</p>

  <ul><li>Pant, R., Russell, T., Glasgow, G., Verschuur, J., Gavin, H., Fowler,
  T. &amp; Hall, J.W. (2021). <em>Analytics for Financial Risk Management of
  Critical Infrastructure in South East Asia – Final Report.</em> Oxford
  Infrastructure Analytics Ltd., Oxford, UK. (Available on request from the
  World Bank)</li></ul>

  <p>The tool being used to visualize the model outputs is developed and
  documented here: </p>

  <ul><li><a href="https://github.com/oi-analytics/seasia"
  target="blank">github.com/oi-analytics/seasia</a></li></ul>

  <p>The outputs specific to Vietnam visualized here were generated from a model
  created and documented here: </p>

  <ul><li><a href="https://github.com/oi-analytics/vietnam-transport"
  target="blank">Vietnam Transport Study</a></li></ul>

  <h2>Funding support</h2>

  <p>This project is led by the Disaster Risk Financing and Insurance Program
  (DRFIP) of the World Bank with support from the Japan&mdash;World Bank Program
  for Mainstreaming DRM in Developing Countries, which is financed by the
  Government of Japan and managed by the Global Facility for Disaster Reduction
  and Recovery (GFDRR) through the Tokyo Disaster Risk Management Hub. </p>

  </article>
);

export default PageIntro;
