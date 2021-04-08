import React from 'react';

const PageIntro = () => (
  <article>
  <h1 className="h1">Southeast Asia Infrastructure Risk: Prototype</h1>

  <p>
    This protoype tool presents data from the World Bank's Transport Risk Study
    of Vietnam within the risk visualisation tool developed as part of the World Bank's
    Argentia Transport Risk Study. Both studies have been undertaken for the World Bank by
    Oxford Infrastructure Analytics Ltd.
  </p>


  <p>
    The purpose of the prototype is to illustrate the type and nature of simulation results and network
    data currently available in the SEA region, and how - using the existing toolchain - this data might
    be accessed and interrogated at a National or Super-National scale.
  </p>


  <p>
    The modelling and analysis presented here aim to support decision-making by identifying
    spatial criticailities, risks, and the performance of adaptation options under current and
    future fluvial flooding outlooks. It comprises a network flow model, generation of failure
    scenarios, economic impact assessment, and cost-benefit analysis of adaptation options.
  </p>

  <p>
    The concepts and model results presented here are documented in the study report:
  </p>

  <p>
    Pant, R., Koks, E.E., Paltan, H., Russell, T., &amp; Hall, J.W. (2019). Argentina – Transport risk analysis.
    Final Report, Oxford Infrastructure Analytics Ltd., Oxford, UK. (Available by request from World Bank)
  </p>


  <p>
    The tool being used to visualize the model outputs was created and documented
    here:
  </p>
  <p>
    <a href="https://github.com/oi-analytics/argentina-transport" target="blank">Argentina Transport Study</a>
  </p>

  <p>
    <a href="https://argentina-transport-risk-analysis.readthedocs.io/en/latest/?badge=latest" target="blank">ReadTheDocs resources</a>
  </p>


  <p>
    The outputs visualized here were generated from a model created and documented
    here:
  </p>

  <p>
    <a href="https://github.com/oi-analytics/vietnam-transport" target="blank">Vietnam Transport Study</a><br></br>
  </p>


  <h1 className="h1">Funding support</h1>
  <p>

  This results inquirer tool has been developed for the Government of Argentina with
  funding support from the World Bank Group and Global Facility for Disaster
  Reduction and Recovery (GFDRR).

  </p>
  </article>
);

export default PageIntro;
