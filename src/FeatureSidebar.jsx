import React, { Fragment } from 'react';

import { commas, titleCase } from './helpers';

const FeatureSidebar = (props) => {
  if (!props.feature) {
    return null
  }
  const f = props.feature.properties;

  const hazard_types = ["fluvial"]
  const scenarios = ["baseline", "rcp_4p5", "rcp_8p5"]
  const hazard_vars = ["height", "probability", "exposure_length"]
  const risk_vars = ["ead","eael_per_day"]
  const adapt_vars = ["ini_adap_cost","maintenance_cost","tot_adap_cost"]
  const adapt_vars_perkm = ["ini_adap_cost_perkm","maintenance_cost_perkm","tot_adap_cost_perkm"]

  return (
    <div className="custom-map-control top-right selected-feature">
      <h4 className="h5">Selected Asset</h4>

      <details>
        <summary>Attributes</summary>
        <dl>
          <dt>ID</dt>
          <dd>{f.node_id || f.edge_id}</dd>
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
            (f.road_class)? (
              <Fragment>
                <dt>Road classification</dt>
                <dd>{titleCase(String(f.road_class))}</dd>
              </Fragment>
            ) : null
          }
          {
            (f.road_cond)? (
              <Fragment>
                <dt>Road Condition</dt>
                <dd>{f.road_cond}</dd>
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
            (f.max_tons && f.min_tons)? (
              <Fragment>
                <dt>Freight flows (tons/day)</dt>
                <dd>{commas(f.min_tons.toFixed(0))} – {commas(f.max_tons.toFixed(0))}</dd>
              </Fragment>
            ) : null
          }
          {
            (f.vehicle_co)? (
              <Fragment>
                <dt>Vehicle counts (vehicles/day)</dt>
                <dd>{commas(f.vehicle_co.toFixed(0))}</dd>
              </Fragment>
            ) : null
          }
        </dl>
      </details>
      {
        (f.flooding_baseline_min_height || f.flooding_rcp_4p5_min_height || f.flooding_rcp_8p5_min_height ||
           f.flooding_baseline_max_height || f.flooding_rcp_4p5_max_height || f.flooding_rcp_8p5_max_height)?
          <details>
            <summary>Flood exposure statistics</summary>
            <dl>
              {
                <table className="table">
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
                                  f["flooding_" + scenario + "_min_" + hazard_var]
                                  && f["flooding_" + scenario + "_max_" + hazard_var]
                                ){
                                  if (hazard_var === "probability") {
                                    if(f["flooding_" + scenario + "_max_" + hazard_var] != null){
                                      return (<td key={hazard_var}>{
                                        `1/${(1 / f["flooding_" + scenario + "_max_" + hazard_var]).toFixed(0)}`
                                      }</td>)
                                    }else{
                                      return(<td key={hazard_var}>-</td>)
                                    }
                                  } else {
                                    return (<td key={hazard_var}>{
                                      commas(f["flooding_" + scenario + "_min_" + hazard_var].toFixed(1))
                                    } – {
                                      commas(f["flooding_" + scenario + "_max_" + hazard_var].toFixed(1))
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
              <dt>No flooding estimated</dt>
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
        :
          <details>
            <summary>Criticality metrics</summary>
            <dl>
              <dt>No values estimated</dt>
            </dl>
          </details>
      }
      {
        (f.baseline_max_ead || f.rcp_4p5_max_ead || f.rcp_8p5_max_ead || 
          f.baseline_max_eael_per_day || f.rcp_4p5_max_eael_per_day || f.rcp_8p5_max_eael_per_day)?
          <details>
            <summary>Risk estimates</summary>
            <dl>
              {
                <table className="table">
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
                              const min_risk = f[scenario + "_min_" + risk_var];
                              const max_risk = f[scenario + "_max_" + risk_var];

                              return (
                                max_risk > min_risk ? 
                                    <td key={risk_var}>{commas(min_risk.toFixed(0))} - {commas(max_risk.toFixed(0))}</td> :
                                    <td key={risk_var}>{commas(min_risk.toFixed(0))}</td>);
                              }
                            )
                          }
                        </tr>)
                      })
                    }
                  </tbody>
                </table>
              }
            </dl>
          </details>
        :
          <details>
            <summary>Risk estimates</summary>
            <dl>
              <dt>No values estimated</dt>
            </dl>
          </details>
      }
      {
        (f.options)?
          <details>
            <summary>Adaptation option</summary>
            <dl>
            {
              <Fragment>
                  <dd>{f.options}</dd>
              </Fragment>
            }
            </dl>
          </details>
          :
          <details>
            <summary>Adaptation option</summary>
            <dl>
              <dt>No values estimated</dt>
            </dl>
          </details>
      }
      {
        (f.baseline_max_ini_adap_cost || f.rcp_4p5_max_ini_adap_cost || f.rcp_8p5_max_ini_adap_cost)?
          <details>
            <summary>Adaptation cost estimates</summary>
            <dl>
              {
                <table className="table">
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
                              if (f[scenario + "_min_" + adapt_var] || f[scenario + "_max_" + adapt_var]) {
                                const min_val = f[scenario + "_min_" + adapt_var];
                                const max_val = f[scenario + "_max_" + adapt_var];
  
                                return (
                                  Math.floor(Math.max(min_val,max_val)) > Math.floor(Math.min(min_val,max_val)) ? 
                                      <td key={adapt_var}>{commas(Math.min(min_val,max_val).toFixed(0))} - {commas(Math.max(min_val,max_val).toFixed(0))}</td> :
                                      <td key={adapt_var}>{commas(min_val.toFixed(0))}</td>
                                );
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
        (f.baseline_max_ini_adap_cost_perkm || f.rcp_4p5_max_ini_adap_cost_perkm || f.rcp_8p5_max_ini_adap_cost_perkm)?
          <details>
            <summary>Adaptation cost estimates per km</summary>
            <dl>
              {
                <table className="table">
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
                              if (f[scenario + "_min_" + adapt_var_perkm] || f[scenario + "_max_" + adapt_var_perkm]) {
                                const min_val = f[scenario + "_min_" + adapt_var_perkm];
                                const max_val = f[scenario + "_max_" + adapt_var_perkm];

                                return (
                                  Math.floor(Math.max(min_val,max_val)) > Math.floor(Math.min(min_val,max_val)) ? 
                                      <td key={adapt_var_perkm}>{commas(Math.min(min_val,max_val).toFixed(0))} - {commas(Math.max(min_val,max_val).toFixed(0))}</td> :
                                      <td key={adapt_var_perkm}>{commas(min_val.toFixed(0))}</td>
                                );
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
        (f.baseline_max_tot_adap_cost || f.rcp_4p5_max_tot_adap_cost || f.rcp_8p5_max_tot_adap_cost)?
          <details>
            <summary>Benefit-cost ratio estimates</summary>
            <BCRWidget
              feature={f}
              scenarios={scenarios}
              updateBCR={props.updateBCR}
              duration={props.duration}
              growth_rate_percentage={props.growth_rate_percentage}
              />
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
      duration: props.duration || 30,
      growth_rate_percentage: props.growth_rate_percentage || 2.8
    }
    this.handleChange = this.handleChange.bind(this);
  }

  handleChange(e) {
    const growth_rate_percentage = (e.target.name === 'growth_rate_percentage')?
      +e.target.value
      : this.state.growth_rate_percentage;

    const duration = (e.target.name === 'duration')?
      +e.target.value
      : this.state.duration;

    const { discount_rate_norm, discount_rate_growth } = discount_rates(growth_rate_percentage / 100)

    this.props.updateBCR({
      duration: duration,
      growth_rate_percentage: growth_rate_percentage,
      discount_growth: discount_rate_growth,
      discount_norm: discount_rate_norm
    })

    this.setState({
      [e.target.name]: +e.target.value  // + coerces to number
    })
  }

  render() {
    const f = this.props.feature;
    const duration = this.state.duration;
    const growth_rate = this.state.growth_rate_percentage / 100;

    return <form>
      <div className="form-group">
        <label htmlFor="duration">Duration of disruption (days)</label>
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
        <label htmlFor="growth_rate_percentage">Growth Rate (%)</label>
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
            const min_ead = f[`${scenario}_min_ead`];
            const max_ead = f[`${scenario}_max_ead`];
            const min_eael_per_day = f[`${scenario}_min_eael_per_day`];
            const max_eael_per_day = f[`${scenario}_max_eael_per_day`];
            const tot_adap_cost = f[`${scenario}_max_tot_adap_cost`];

            if (!min_ead || !max_ead || !min_eael_per_day || !max_eael_per_day || !tot_adap_cost) {
              return null
            }

            const data = calculateAdaption(
              (max_ead + min_ead) / 2.0, min_eael_per_day, max_eael_per_day, tot_adap_cost, duration, growth_rate)

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
  const { discount_rate_norm, discount_rate_growth } = discount_rates(growth_rate)
  const min_benefit = benefit(discount_rate_norm, ead, duration, discount_rate_growth, min_eael_per_day)
  const max_benefit = benefit(discount_rate_norm, ead, duration, discount_rate_growth, max_eael_per_day)
  const min_bcr = bcr(min_benefit, tot_adap_cost)
  const max_bcr = bcr(max_benefit, tot_adap_cost)

  return { min_benefit, min_bcr, max_benefit, max_bcr }
}

function discount_rates(growth_rate) {
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

  return { discount_rate_norm, discount_rate_growth }
}

function benefit(discount_rate_norm, ead, duration, discount_rate_growth, eael_per_day){
  return (discount_rate_norm * ead) + (duration * discount_rate_growth * eael_per_day)
}

function bcr(benefit, tot_adap_cost){
  return benefit / tot_adap_cost;
}

export default FeatureSidebar;
