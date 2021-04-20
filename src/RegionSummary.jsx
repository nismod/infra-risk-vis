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
      <h2>Annual direct damages and economic losses, historical scenario</h2>
      <VegaLite
        spec={{
          width: 400,
          height: 200,
          mark: 'area',
          encoding: {
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
          data: { url: 'historical_total/historical_total_'+r.GID_1+'.csv' },
        }}
        // defines the actions available behind the ... menu top-right of chart
        actions={{
          export: true,
          source: false,
          compiled: false,
          editor: false,
        }}
      />
      <h2>Max direct damages across climate scenarios (historical/2080)</h2>
      <VegaLite
        spec={{
          width: 400,
          height: 200,
          mark: 'line',
          encoding: {
            x: { field: 'probability', type: 'quantitative', title: 'Probability'},
            y: { field: 'assetDamage', type: 'quantitative', title: 'Annual Damages (US$m)'},
            color: { field: 'rcp', type: 'nominal', title: 'RCP'},
            row: { field: 'hazard'}
          },
          data: { url: 'climate_total/climate_total_'+r.GID_1+'.csv' },
        }}
        // defines the actions available behind the ... menu top-right of chart
        actions={{
          export: true,
          source: false,
          compiled: false,
          editor: false,
        }}
      />
      <h2>Max indirect damages across climate scenarios</h2>
      <VegaLite
        spec={{
          width: 400,
          height: 200,
          mark: 'line',
          encoding: {
            x: { field: 'probability', type: 'quantitative', title: 'Probability'},
            y: { field: 'gdp', type: 'quantitative', title: 'Annual Damages (US$m)'},
            color: { field: 'rcp', type: 'nominal', title: 'RCP'},
            row: { field: 'hazard'}
          },
          data: { url: 'climate_total/climate_total_'+r.GID_1+'.csv' },
        }}
        // defines the actions available behind the ... menu top-right of chart
        actions={{
          export: true,
          source: false,
          compiled: false,
          editor: false,
        }}
      />
      <h2>Total expected damages and losses under climate scenarios</h2>
      <p>* Note that cyclone damages are only estimated for historical scenario.</p>
      <table className="table">
        <thead>
          <tr>
            <th>Climate Scenario</th>
            <th>Expected Annual Damages (US$m)</th>
            <th>Expected Annual Economic Losses (30-day disruption, US$m)</th>
          </tr>
        </thead>
        <tbody>
          <tr>
            <th>Historical</th>
            <td>{r.minEAD.toFixed(1)}&ndash;{r.maxEAD.toFixed(1)}</td>
            <td>{(r['EAEL-gdp'] * disruption_fraction).toFixed(1)}</td>
          </tr>
          <tr>
            <th>RCP 4.5</th>
            <td>{r.minEAD_rcp4p5.toFixed(1)}&ndash;{r.maxEAD_rcp4p5.toFixed(1)}</td>
            <td>{(r['EAEL-gdp_rcp4p5'] * disruption_fraction).toFixed(1)}</td>
          </tr>
          <tr>
            <th>RCP 8.5</th>
            <td>{r.minEAD_rcp8p5.toFixed(1)}&ndash;{r.maxEAD_rcp8p5.toFixed(1)}</td>
            <td>{(r['EAEL-gdp_rcp8p5'] * disruption_fraction).toFixed(1)}</td>
          </tr>
        </tbody>
      </table>
    </article>
  );
}

export default RegionSummary;
