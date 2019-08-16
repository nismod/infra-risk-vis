import React, { Fragment } from 'react';

import { commas, titleCase } from './helpers';

const FeatureSidebar = (props) => {
  if (!props.feature) {
    return null
  }
  const f = props.feature.properties;

  const hazard_types = ["pluvial_flooding", "fluvial_flooding"]
  const scenarios = ["baseline", "future_med", "future_high"]
  const hazard_vars = ["flood_depth", "probability", "exposure_length"]
  const risk_vars = ["ead","eael_per_day"]
  const adapt_vars = ["ini_adap_cost","tot_maintenance_cost","tot_adap_cost"]
  const adapt_vars_perkm = ["ini_adap_cost_per_km","tot_maintenance_cost_per_km","tot_adap_cost_per_km"]

  return (
    <div className="custom-map-control top-right selected-feature">
      <h4 className="h5">Selected Asset</h4>

      <details>
        <summary>Attributes</summary>
        <dl>
          <dt>ID</dt>
          <dd>{f.node_id || f.edge_id || f.bridge_id}</dd>
          {
            (f.name)? (
              <Fragment>
                <dt>Name</dt>
                <dd>{f.name.replace("0","Unknown")}</dd>
              </Fragment>
            ) : null
          }
          {
            (f.road_name)? (
              <Fragment>
                <dt>Road name</dt>
                <dd>{f.road_name.replace(",","/")}</dd>
              </Fragment>
            ) : null
          }
          {
            (f.ruta)? (
              <Fragment>
                <dt>Road name</dt>
                <dd>{f.ruta}</dd>
              </Fragment>
            ) : null
          }
          {
            (f.iata)? (
              <Fragment>
                <dt>IATA code</dt>
                <dd>{f.iata}</dd>
              </Fragment>
            ) : null
          }
          {
            (f.locality)? (
              <Fragment>
                <dt>Port cluster</dt>
                <dd>{f.locality}</dd>
              </Fragment>
            ) : null
          }
          {
            (f.operador)? (
              <Fragment>
                <dt>Line Operating company</dt>
                <dd>{f.operador}</dd>
              </Fragment>
            ) : null
          }
          {
            (f.linea)? (
              <Fragment>
                <dt>Line name</dt>
                <dd>{f.linea}</dd>
              </Fragment>
            ) : null
          }
          {
            (f.road_type)? (
              <Fragment>
                <dt>Road classification</dt>
                <dd>{titleCase(f.road_type)}</dd>
              </Fragment>
            ) : null
          }
          {
            (f.structure_type)? (
              <Fragment>
                <dt>Structure type</dt>
                <dd>{f.structure_type}</dd>
              </Fragment>
            ) : null
          }
          {
            (f.pavement_material_asc)? (
              <Fragment>
                <dt>Existing pavement material</dt>
                <dd>{f.pavement_material_asc.replace("0","Unknown").replace(",","/")}</dd>
              </Fragment>
            ) : null
          }
          {
            (f.surface)? (
              <Fragment>
                <dt>Existing pavement material</dt>
                <dd>{f.surface.replace("0","Unknown").replace(",","/")}</dd>
              </Fragment>
            ) : null
          }
          {
            (f.length)? (
              <Fragment>
                <dt>Link length (km)</dt>
                <dd>{f.length.toFixed(2)}</dd>
              </Fragment>
            ) : null
          }
          {
            (f.width)? (
              <Fragment>
                <dt>Width (m)</dt>
                <dd>{f.width.toFixed(2)}</dd>
              </Fragment>
            ) : null
          }
          {
            (f.min_speed && f.max_speed)? (
              <Fragment>
                <dt>Speeds (km/hr)</dt>
                <dd>{f.min_speed.toFixed(0)} – {f.max_speed.toFixed(0)}</dd>
              </Fragment>
            ) : null
          }
          {
            (f.max_total_tons && f.min_total_tons)? (
              <Fragment>
                <dt>Freight flows (tons/day)</dt>
                <dd>{commas(f.min_total_tons.toFixed(0))} – {commas(f.max_total_tons.toFixed(0))}</dd>
              </Fragment>
            ) : null
          }
          {
            (f.passengers)? (
              <Fragment>
                <dt>Passengers</dt>
                <dd>{commas(f.passengers.toFixed(0))}</dd>
              </Fragment>
            ) : null
          }
          {
            (f.tmda_count)? (
              <Fragment>
                <dt>Vehicle counts (Annual)</dt>
                <dd>{commas(f.tmda_count.toFixed(0))}</dd>
              </Fragment>
            ) : null
          }
        </dl>
      </details>
      {
        (f.pluvial_flooding_baseline_min_flood_depth || f.fluvial_flooding_baseline_min_flood_depth || f.pluvial_flooding_future_med_min_flood_depth || f.fluvial_flooding_future_med_min_flood_depth || f.pluvial_flooding_future_high_min_flood_depth || f.fluvial_flooding_future_high_min_flood_depth)?
          <details>
            <summary>Flood exposure statistics</summary>
            <dl>
              {
                <table class="table">
                  <thead>
                    <tr>
                      <th>Flood type</th>
                      <th>Climate scenario</th>
                      <th>Flood Depth (m)</th>
                      <th>Exceedance probability (1/years)</th>
                      <th>Exposure length (m)</th>
                    </tr>
                  </thead>
                  <tbody>
                    {
                      hazard_types.map((hazard) => {
                        return scenarios.map((scenario) => {
                          return (<tr key={scenario + hazard}>
                            <td>{titleCase(hazard.replace("_flooding", ""))}</td>
                            <td>{titleCase(scenario.replace("_"," "))}</td>
                            {
                              hazard_vars.map((hazard_var) => {
                                if (
                                  f[hazard + "_" + scenario + "_min_" + hazard_var]
                                  && f[hazard + "_" + scenario + "_max_" + hazard_var]
                                ){
                                  if (hazard_var === "probability") {
                                    return (<td key={hazard_var}>{
                                      `1/${(1 / f[hazard + "_" + scenario + "_max_" + hazard_var]).toFixed(0)}`
                                    }</td>)
                                  } else {
                                    return (<td key={hazard_var}>{
                                      commas(f[hazard + "_" + scenario + "_min_" + hazard_var].toFixed(1))
                                    } – {
                                      commas(f[hazard + "_" + scenario + "_max_" + hazard_var].toFixed(1))
                                    }</td>)
                                  }
                                } else {
                                    return (<td key={hazard_var}>-</td>)
                                }
                              })
                            }
                          </tr>)
                        })
                      })
                    }
                  </tbody>
                </table>
              }
            </dl>
          </details>
        :
          <details>
            <summary>Flood exposure statistics</summary>
            <dl>
              <dt>No values</dt>
            </dl>
          </details>
      }
      {
        (f.max_tr_loss || f.max_econ_loss || f.max_econ_impact)?
          <details>
            <summary>Criticality metrics</summary>
            <dl>
              {
                (f.max_tr_loss && f.min_tr_loss)? (
                  <Fragment>
                    <dt>Rerouting loss (USD/day)</dt>
                    <dd>{commas(f.min_tr_loss.toFixed(0))} – {commas(f.max_tr_loss.toFixed(0))}</dd>
                  </Fragment>
                ) : null
              }
              {
                (f.max_econ_loss && f.min_econ_loss)? (
                  <Fragment>
                    <dt>Macroeconomic loss (USD/day)</dt>
                    <dd>{commas(f.min_econ_loss.toFixed(0))} – {commas(f.max_econ_loss.toFixed(0))}</dd>
                  </Fragment>
                ) : null
              }
              {
                (f.max_econ_impact && f.min_econ_impact)? (
                  <Fragment>
                    <dt>Total Economic impact (USD/day)</dt>
                    <dd>{commas(f.min_econ_impact.toFixed(0))} – {commas(f.max_econ_impact.toFixed(0))}</dd>
                  </Fragment>
                ) : null
              }
            </dl>
          </details>
        : null
      }
      {
        (f.baseline_ead || f.future_med_ead || f.future_high_ead || f.baseline_max_eael_per_day || f.future_med_max_eael_per_day || f.future_high_max_eael_per_day)?
          <details>
            <summary>Risk estimates</summary>
            <dl>
              {
                <table class="table">
                  <thead>
                    <tr>
                      <th>Climate scenario</th>
                      <th>Expected Annual Damages (USD)</th>
                      <th>Expected Annual Economic Losses/day (USD/day)</th>
                    </tr>
                  </thead>
                  <tbody>
                    {
                      scenarios.map((scenario) => {
                        return (<tr key={scenario}>
                          <td>{titleCase(scenario.replace("_"," "))}</td>
                          {
                            risk_vars.map((risk_var) => {
                              if (risk_var === "ead" && f[scenario + "_" + risk_var]) {
                                return (<td key={risk_var}>{
                                  commas(f[scenario + "_" + risk_var].toFixed(0))
                                }</td>)
                              } else if (f[scenario + "_min_" + risk_var] && f[scenario + "_max_" + risk_var]) {

                                return (<td key={risk_var}>{
                                  commas(f[scenario + "_min_" + risk_var].toFixed(0))
                                } – {
                                  commas(f[scenario + "_max_" + risk_var].toFixed(0))
                                }</td>)
                              } else {
                                return (<td key={risk_var}>-</td>)
                              }
                            })
                          }
                        </tr>)
                      })
                    }
                  </tbody>
                </table>
              }
            </dl>
          </details>
        : null
      }
      {
        (f.baseline_ini_adap_cost || f.future_med_ini_adap_cost || f.future_high_ini_adap_cost)?
          <details>
            <summary>Adaptation cost estimates</summary>
            <dl>
              {
                <table class="table">
                  <thead>
                    <tr>
                      <th>Climate scenario</th>
                      <th>Initial investment (USD)</th>
                      <th>Maintenance investments (USD)</th>
                      <th>Total investments (USD)</th>
                    </tr>
                  </thead>
                  <tbody>
                    {
                      scenarios.map((scenario) => {
                        return (<tr key={scenario}>
                          <td>{titleCase(scenario.replace("_"," "))}</td>
                          {
                            adapt_vars.map((adapt_var) => {
                              if (f[scenario + "_" + adapt_var]) {
                                return (<td key={adapt_var}>{
                                  commas(f[scenario + "_" + adapt_var].toFixed(0))
                                }</td>)
                              } else {
                                return (<td key={adapt_var}>-</td>)
                              }
                            })
                          }
                        </tr>)
                      })
                    }
                  </tbody>
                </table>
              }
            </dl>
          </details>
        : null
      }
      {
        (f.baseline_ini_adap_cost || f.future_med_ini_adap_cost || f.future_high_ini_adap_cost)?
          <details>
            <summary>Adaptation cost estimates per km</summary>
            <dl>
              {
                <table class="table">
                  <thead>
                    <tr>
                      <th>Climate scenario</th>
                      <th>Initial investment (USD)</th>
                      <th>Maintenance investments (USD)</th>
                      <th>Total investments (USD)</th>
                    </tr>
                  </thead>
                  <tbody>
                    {
                      scenarios.map((scenario) => {
                        return (<tr key={scenario}>
                          <td>{titleCase(scenario.replace("_"," "))}</td>
                          {
                            adapt_vars_perkm.map((adapt_var_perkm) => {
                              if (f[scenario + "_" + adapt_var_perkm]) {
                                return (<td key={adapt_var_perkm}>{
                                  commas(f[scenario + "_" + adapt_var_perkm].toFixed(0))
                                }</td>)
                              } else {
                                return (<td key={adapt_var_perkm}>-</td>)
                              }
                            })
                          }
                        </tr>)
                      })
                    }
                  </tbody>
                </table>
              }
            </dl>
          </details>
        : null
      }
      {
        (f.baseline_ini_adap_cost || f.future_med_ini_adap_cost || f.future_high_ini_adap_cost)?
          <details>
            <summary>Benefit-cost ratio estimates</summary>
            <BCRWidget feature={f} scenarios={scenarios} />
          </details>
          : null
      }
    </div>
  )
}

class BCRWidget extends React.Component {
  constructor(props) {
    super(props)
    this.state = {
      duration: 10,
      growth_rate_percentage: 1
    }
    this.handleChange = this.handleChange.bind(this);
  }

  handleChange(e) {
    this.setState({
      [e.target.name]: e.target.value
    })
  }

  render() {
    const f = this.props.feature;
    const duration = this.state.duration;
    const growth_rate = this.state.growth_rate_percentage / 100;

    return <form>
      <div className="form-group">
        <label for="duration">Duration of disruption (days)</label>
        <input
          id="duration"
          name="duration"
          onChange={this.handleChange}
          className="form-control"
          type="range"
          value={this.state.duration}
          min={10}
          max={100}
          step={10}
          list="duration_ticks"
          />
        <datalist id="duration_ticks">
          <option value="10" label="10" />
          <option value="20" label="20" />
          <option value="30" label="30" />
          <option value="40" label="40" />
          <option value="50" label="50" />
          <option value="60" label="60" />
          <option value="70" label="70" />
          <option value="80" label="80" />
          <option value="90" label="90" />
          <option value="100" label="100" />
        </datalist>
        <span>{this.state.duration}</span>
      </div>
      <div className="form-group">
        <label for="growth_rate_percentage">Growth Rate (%)</label>
        <input
          id="growth_rate_percentage"
          name="growth_rate_percentage"
          onChange={this.handleChange}
          className="form-control"
          type="range"
          value={this.state.growth_rate_percentage}
          min={-2}
          max={4}
          step={0.2}
          list="growth_rate_percentage_ticks"
          />
        <datalist id="growth_rate_percentage_ticks">
          <option value="-2" label="-2%" />
          <option value="-1" />
          <option value="0" label="0%" />
          <option value="1" />
          <option value="2" label="2%" />
          <option value="3" />
          <option value="4" label="4%" />
        </datalist>
        <span>{this.state.growth_rate_percentage.toFixed(1)}</span>
      </div>
      <table className="table">
        <thead>
          <tr>
            <th>Scenario</th>
            <th>Cost</th>
            <th>Benefit</th>
            <th>Benefit-Cost Ratio</th>
          </tr>
        </thead>
        <tbody>
        {
          this.props.scenarios.map(scenario => {
            // data from asset
            const ead = f[`${scenario}_ead`];
            const min_eael_per_day = f[`${scenario}_min_eael_per_day`];
            const max_eael_per_day = f[`${scenario}_max_eael_per_day`];
            const tot_adap_cost = f[`${scenario}_tot_adap_cost`];

            if (!ead || !min_eael_per_day || !max_eael_per_day || !tot_adap_cost) {
              return null
            }

            const data = calculateAdaption(
              ead, min_eael_per_day, max_eael_per_day, tot_adap_cost, duration, growth_rate)

            return <tr key={scenario}>
              <td>{titleCase(scenario.replace('_', ' '))}</td>
              <td>{commas(tot_adap_cost.toFixed(0))}</td>
              <td>{`${commas(data.min_benefit.toFixed(0))} – ${commas(data.max_benefit.toFixed(0))}`}</td>
              <td>{`${data.min_bcr.toFixed(2)} – ${data.max_bcr.toFixed(2)}`}</td>
            </tr>
          })
        }
        </tbody>
      </table>
    </form>
  }
}

/**
 * Calculate min/max adaptation benefits (as avoided damages and macroeconomic losses) and bcr,
 * given hazard duration and economic growth rate assumptions.
 *
 * @param {number} ead from asset - expected annual damages
 * @param {number} min_eael_per_day from asset - min expected annual losses per day
 * @param {number} max_eael_per_day from asset - max expected annual losses per day
 * @param {number} tot_adap_cost from asset - total adaptation cost
 * @param {number} duration as integer number of days
 * @param {number} growth_rate as rate between 0 and 1 (typically presented as percentage)
 */
function calculateAdaption(ead, min_eael_per_day, max_eael_per_day, tot_adap_cost, duration, growth_rate) {
  // See atra.adaptation_options.calculate_discounting_arrays() at
  // https://github.com/oi-analytics/argentina-transport/blob/master/src/atra/adaptation_options.py
  // for reference.

  // fix parameters
  const discount_rate = 0.12;
  const start_year = 2016;
  const end_year = 2050;

  const years_from_start = Array.from(
    Array(end_year - start_year).keys()
  )

  const discount_rate_norms = years_from_start.map(
    year => 1.0 / Math.pow(1.0 + discount_rate, year)
  )
  const discount_rate_norm = discount_rate_norms.reduce((a, b) => a + b, 0)

  const discount_rate_growths = years_from_start.map(
    year => Math.pow(1.0 + growth_rate, year) /
            Math.pow(1.0 + discount_rate, year)
  )
  const discount_rate_growth = discount_rate_growths.reduce((a, b) => a + b, 0)

  console.log(years_from_start)
  console.log(discount_rate_norms, discount_rate_norm)
  console.log(discount_rate_growths, discount_rate_growth)

  // for variable growth_rate and duration we want to estimate min/max benefit and bcr
  const min_benefit = benefit(discount_rate_norm, ead, duration, discount_rate_growth, min_eael_per_day)
  const max_benefit = benefit(discount_rate_norm, ead, duration, discount_rate_growth, max_eael_per_day)
  const min_bcr = bcr(min_benefit, tot_adap_cost)
  const max_bcr = bcr(max_benefit, tot_adap_cost)

  return { min_benefit, min_bcr, max_benefit, max_bcr }
}

function benefit(discount_rate_norm, ead, duration, discount_rate_growth, eael_per_day){
  return (discount_rate_norm * ead) + (duration * discount_rate_growth * eael_per_day)
}

function bcr(benefit, tot_adap_cost){
  return benefit / tot_adap_cost;
}

export default FeatureSidebar;
