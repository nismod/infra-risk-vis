import React from 'react';

const PageIntro = () => (
  <article>
  <h1>Jamaica Infrastructure Risk Assessment Tool</h1>

  <p>This prototype tool presents infrastructure risk analytics for Jamaica.</p>

  <p>The modelling and analysis presented here aim to support climate adaptation
  decision-making by identifying spatial criticalities and risks under current
  and future hazard scenarios. It comprises a direct damage estimation and an
  indirect economic loss estimation of GDP disruptions due to asset failures and
  service disruption.</p>

  <table className="table table-sm table-striped">
    <thead>
      <tr>
        <th>Infrastructure</th>
        <th>Assets</th>
        <th>Expected Annual Damages (EAD)</th>
        <th>Expected Annual Economic Losses (EAEL)</th>
      </tr>
    </thead>
    <tbody>
      <tr>
        <td rowspan="2">Transport</td>
        <td>Road links</td>
        <td>Cost of rehabilitation/reinstating damaged assets</td>
        <td>Rerouting costs + wider effects of service disruption</td>
      </tr>
      <tr>
        <td>Railway tracks</td>
        <td>Cost of rehabilitation/reinstating damaged assets</td>
        <td>Rerouting costs + wider effects of service disruption</td>
      </tr>
      <tr>
        <td>Energy</td>
        <td>Electricity transmission and distribution grid: generation, lines and substations</td>
        <td>Cost of rehabilitation/reinstating damaged assets</td>
        <td>Wider effects of service disruption</td>
      </tr>
      <tr>
        <td>Water</td>
        <td>Water supply and wastewater networks, abstraction points, irrigation schemes</td>
        <td>Cost of rehabilitation/reinstating damaged assets</td>
        <td>Wider effects of service disruption</td>
      </tr>
    </tbody>
  </table>

  <p>This tool to visualize the model outputs is developed and documented
  here:</p>

  <ul><li><a href="https://github.com/nismod/infra-risk-vis"
  target="blank">github.com/nismod/infra-risk-vis</a></li></ul>

  <p>The analytics for Jamaica are produced using the code and models here:</p>

  <ul><li><a href="https://github.com/nismod/jamaica-infrastructure"
  target="blank">github.com/nismod/jamaica-infrastructure</a></li></ul>

  <h2>Funding support</h2>

  <p>This research is led by researchers in the Oxford Programme for Sustainable
  Infrastructure Systems in the Environmental Change Institute, University of
  Oxford, for the Government of Jamaica (GoJ) as part of a project funded by UK
  Aid (FCDO). The initiative forms part of the Coalition for Climate Resilient
  Investment’s (CCRI) collaboration with the GoJ, which includes analysis of
  nature-based approaches to build resilience in Jamaica to be procured and
  funded by the Green Climate Fund (GCF).</p>

  </article>
);

export default PageIntro;
