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
                <dd>{f.min_speed.toFixed(0)}-{f.max_speed.toFixed(0)}</dd>
              </Fragment>
            ) : null
          }
          {
            (f.max_total_tons && f.min_total_tons)? (
              <Fragment>
                <dt>Freight flows (tons/day)</dt>
                <dd>{commas(f.min_total_tons.toFixed(0))}-{commas(f.max_total_tons.toFixed(0))}</dd>
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
                                    }-{
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
                    <dd>{commas(f.min_tr_loss.toFixed(0))}-{commas(f.max_tr_loss.toFixed(0))}</dd>
                  </Fragment>
                ) : null
              }
              {
                (f.max_econ_loss && f.min_econ_loss)? (
                  <Fragment>
                    <dt>Macroeconomic loss (USD/day)</dt>
                    <dd>{commas(f.min_econ_loss.toFixed(0))}-{commas(f.max_econ_loss.toFixed(0))}</dd>
                  </Fragment>
                ) : null
              }
              {
                (f.max_econ_impact && f.min_econ_impact)? (
                  <Fragment>
                    <dt>Total Economic impact (USD/day)</dt>
                    <dd>{commas(f.min_econ_impact.toFixed(0))}-{commas(f.max_econ_impact.toFixed(0))}</dd>
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
                                }-{
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
    </div>
  )
}

export default FeatureSidebar;
