import React from 'react';

// TODO: Change the words - Raghav to advise?
const FloodHelp = () => (
  <div className="custom-map-control top-right selected-feature">
    <h4 className="h5">Flood Climate Outlooks - Explanation</h4>
    <dl>
      <dt>Baseline</dt>
      <dd>{"The estimated flooded depths and areas averaged over 1960-1999 based on historical rainfall records."}</dd>
      <dt>Future Median</dt>
      <dd>{" The estimated flooded depths and areas averaged over 2010-2049 based on global climate model outputs assuming global carbon emissions peak by 2040 before declining. The layer displayed here shows the flood outlines from the UK Met Office Hadley Centre Global Environment Model version 2 (HadGEM2-ES) model output, which is 1 of 5 models used in this study."}</dd>
      <dt>Future High</dt>
      <dd>{" The estimated flooded depths and areas averaged over 2010-2049 based on global climate model outputs assuming global carbon emissions continue to rise throughout the 21st century. The layer displayed here shows the flood outlines from the UK Met Office Hadley Centre Global Environment Model version 2 (HadGEM2-ES) model output, which is 1 of 5 models used in this study."}</dd>
    </dl>
  </div>
)

export default FloodHelp;
