import React from 'react';

const FloodHelp = () => (
  <div className="custom-map-control top-right selected-feature">
    <h4 className="h5">Flood Climate Outlooks - Explanation</h4>
    <dl>
      <dt>Baseline</dt>
      <dd>{"The estimated flooded depths and areas informed by the 1986-2005 baseline precipitation given by historical rainfall records."}</dd>
      <dt>Future Median</dt>
      <dd>{"The estimated flooded depths and areas by 2050 informed by the median value of 5-day maximum precipitation change with respect to the 1986-2005 baseline precipitation, which was approximately +6%, recorded across 32 GCM's and across RCP 4.5 and RCP 8.5 climate emission scenarios."}</dd>
      <dt>Future High</dt>
      <dd>{"The estimated flooded depths and areas by 2050 informed by 90th percentile value of 5-day maximum precipitation changes with respect to the 1986-2005 baseline precipitation, which was approximately +12%, recorded across 32 GCM's  and across RCP 4.5 and RCP 8.5 climate emission scenarios."}</dd>
    </dl>
  </div>
)

export default FloodHelp;
