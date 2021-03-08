import React from 'react'
import { commas, titleCase } from './helpers'

const Tooltip = (props) => {
  const entries = {}

  for (const f of props.features) {
    let title = titleCase(
      f.sourceLayer.replace(/_/g," ")
        .replace(/m(\d)/, 'm-$1')
        .replace('4m-999m', '>4m')
        .replace('1in', '1/')
        .replace('edges', '')
    );
    let subtitle = (f.properties.road_type)?
      "(" + f.properties.road_type + ")"
      : "";

    let max_value;
    let detail;

    if (props.map_style === "roads" || props.map_style === "rail" || props.map_style === "airwater" || props.map_style === "overview") {
      if (f.sourceLayer === "air_nodes") {
        max_value = f.properties.passengers;
        detail = (f.properties.passengers)?
          commas(f.properties.passengers.toFixed(0)) + " passengers"
          : "";

      } else {
        max_value = f.properties.max_tons;

        detail = (f.properties.max_tons && f.properties.min_tons)?
          " " +
          commas(f.properties.min_tons.toFixed(0)) + " – " +
          commas(f.properties.max_tons.toFixed(0)) + " tons/day freight flows"
          : "";
      }
    }

    if (props.map_style === "energy_network") {
      detail = "Asset Ref: " + f.properties.fid;
    }

    if (props.map_style === "impact") {
      max_value = f.properties.max_econ_impact;

      detail = (f.properties.max_econ_impact && f.properties.min_econ_impact)?
        " " +
        commas(f.properties.min_econ_impact.toFixed(0)) + " – " +
        commas(f.properties.max_econ_impact.toFixed(0)) + " USD/day total economic impact"
        : "";
    }

    if (props.map_style === "risk") {
      max_value = Math.max(
        ((f.properties.baseline_max_ead || 0) + (f.properties.baseline_max_eael_per_day || 0) *30),
        ((f.properties.rcp_4p5_max_ead || 0) + (f.properties.rcp_4p5_max_eael_per_day || 0) * 30),
        ((f.properties.rcp_8p5_max_ead || 0) + (f.properties.rcp_8p5_max_eael_per_day || 0) * 30)
      );

      detail = (max_value)?
        " up to " + commas(max_value.toFixed(0)) + " USD expected annual damages plus losses for a 30-day disruption"
        : "";
    }

    if (!entries[f.sourceLayer] || entries[f.sourceLayer].max_value < max_value) {
      entries[f.sourceLayer] = { title, subtitle, max_value, detail }
    }
  }

  return (props.features.length)? (
    <div className="tooltip-wrap">
      <div className="tooltip-body">
        {
          Object.values(entries).map((entry, i) => {
            return (
            <div key={i}>
              <strong>{entry.title} {entry.subtitle}</strong>
              { entry.detail }
            </div>
            )
          })
        }
      </div>
      <span className="tooltip-triangle"></span>
    </div>
  ) : null;
}

export default Tooltip;
