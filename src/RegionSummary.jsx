import React from 'react';
import  { VegaLite } from 'react-vega';

const RegionSummary = (props) => {
  if (!props.region) {
    return (
      <article>
        <h1 className="h2">Select a region</h1>
        <p className="alert alert-primary">
          Click the map to see details for a region.
        </p>
      </article>
    );
  }
  console.log(props.region);
  return (
    <article>
      <h1 className="h2">{props.region.ISO_A3 || props.region.ISO3_CODE}</h1>
      <p>Annual damages and economic losses</p>
      <VegaLite
        spec={{
          width: 400,
          height: 200,
          mark: 'area',
          encoding: {
            // probability,min_econ_impact,max_econ_impact,damages,total_min,total_max,ini_investment
            x: { field: 'probability', type: 'quantitative', title: 'Probability' },
            y: { field: 'value', "aggregate": "sum", title: 'Value ($USDm)' },
            color: { field: 'field', type: 'nominal', title: 'Value'}
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
