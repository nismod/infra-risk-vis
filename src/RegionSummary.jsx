import React from 'react';
import  { VegaLite } from 'react-vega';

const RegionSummary = (props) => {
  // if (!props.region) {
  //   return (
  //     <article>
  //       <h1 className="h2">Select a region</h1>
  //       <p className="alert alert-primary">
  //         Click the map to see details for a region.
  //       </p>
  //     </article>
  //   );
  // }
  return (
    <article>
      {
        (props.region && props.region.NAME_1 && props.region.NAME_0)?
          <p className="alert alert-info">
            {props.region.NAME_1}, {props.region.NAME_0}
          </p>
        : null
      }

      <h1 className="h2">Vietnam</h1>

      <p className="alert alert-primary">Coming soon: admin-1 regional summary
      statistics to be loaded when the map is clicked. Currently showing
      summary statistics as calculated for Vietnam.</p>

      <p>Annual damages and economic losses</p>
      <VegaLite
        spec={{
          width: 400,
          height: 200,
          mark: 'area',
          encoding: {
            // probability,min_econ_impact,max_econ_impact,damages,total_min,total_max,ini_investment
            x: { field: 'probability', type: 'quantitative', title: 'Probability' },
            y: {
              field: 'value',
              aggregate: "sum",
              type: 'quantitative',
              title: 'Value ($USDm)'
            },
            color: {
              field: 'field',
              type: 'nominal',
              title: 'Value',
              "scale": {
                "domain": ["1. Economic impact", "2. Direct damages"],
                "range": ["#1f77b4", "#e7ba52"]
              },
            }
          },
          data: { url: 'aggregated_stats_national_summary.csv' },
        }}
        // defines the actions available behind the ... menu top-right of chart
        actions={{
          export: true,
          source: false,
          compiled: false,
          editor: false,
        }}
        />
      <p>Annual damages across climate scenarios</p>
      <VegaLite
        spec={{
          width: 400,
          height: 200,
          mark: 'line',
          encoding: {
            facet: {
              field: "climate_scenario",
              title: "Climate Scenario",
              type: "nominal",
              columns: 1
            },
            // probability,min_econ_impact,max_econ_impact,damages,total_min,total_max,ini_investment
            x: { field: 'probability', type: 'quantitative', title: 'Probability'},
            y: { field: 'damages', type: 'quantitative', title: 'Annual Damages ($USDm)'},
            color: { field: 'model', type: 'nominal', title: 'Global Climate Model'}
          },
          data: { url: 'aggregated_stats_national.csv' },
        }}
        // defines the actions available behind the ... menu top-right of chart
        actions={{
          export: true,
          source: false,
          compiled: false,
          editor: false,
        }}
        />
    </article>
  );
}

export default RegionSummary;
