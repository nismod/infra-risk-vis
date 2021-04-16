import React from 'react';
import  { VegaLite } from 'react-vega';

const RegionSummary = (props) => {
  if (!props.region) {
    return (
      <article>
        <h1>Region, Country</h1>
        <small>Region code: &hellip;</small>
        <p className="alert alert-primary">
          Click the map to see details for a region.
        </p>
      </article>
    );
  }
  const r = props.region;
  const disruption_fraction = 30 / 365;
  return (
    <article>
      <h1>{r.NAME_1}, {r.NAME_0}</h1>
      <small>Region code: {r.GID_1}</small>
      <h2>Total expected damages and losses under climate scenarios</h2>
      <table className="table">
        <thead>
          <tr>
            <th>Climate Scenario</th>
            <th>Expected Annual Damages</th>
            <th>Expected Annual Economic Losses (30-day disruption)</th>
          </tr>
        </thead>
        <tbody>
          <tr>
            <th>Historical</th>
            <td>{r.minEAD.toFixed(0)}&ndash;{r.maxEAD.toFixed(0)}</td>
            <td>{(r['EAEL-gdp'] * disruption_fraction).toFixed(0)}</td>
          </tr>
          <tr>
            <th>RCP 4.5</th>
            <td>{r.minEAD_rcp4p5.toFixed(0)}&ndash;{r.maxEAD_rcp4p5.toFixed(0)}</td>
            <td>{(r['EAEL-gdp_rcp4p5'] * disruption_fraction).toFixed(0)}</td>
          </tr>
          <tr>
            <th>RCP 8.5</th>
            <td>{r.minEAD_rcp8p5.toFixed(0)}&ndash;{r.maxEAD_rcp8p5.toFixed(0)}</td>
            <td>{(r['EAEL-gdp_rcp8p5'] * disruption_fraction).toFixed(0)}</td>
          </tr>
        </tbody>
      </table>
      {/*
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
              title: 'Value (US$m)'
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
            y: { field: 'damages', type: 'quantitative', title: 'Annual Damages (US$m)'},
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
        */}
    </article>
  );
}

export default RegionSummary;
