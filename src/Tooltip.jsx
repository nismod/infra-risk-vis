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

    if (props.map_style === "roads" || props.map_style === "rail" || props.map_style === "electricity" || props.map_style === "overview") {

      title = "ID: " + (f.properties.osm_id || f.properties.link);

      if (f.properties.EAD_min_usd && f.properties.EAD_max_usd) {
        max_value = f.properties.EAD_max || 0;

        detail = " EAD: $" +
          f.properties.EAD_min_usd + "–" +
          f.properties.EAD_max_usd;

        if (f.properties.EAEL_daily_usd) {
          detail += ", EAEL: $" +
          f.properties.EAEL_daily_usd + ".";
        } else {
          detail += ".";
        }
      } else {
        detail = "(no exposure calculated)  "
      }

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

    // Hazard details
    if (f.properties.depth_m) {
      detail = "Depth: " + f.properties.depth_m.toFixed(1) + "m";
    }
    if (f.properties["gust_speed_ms-1"]) {
      detail = "Gust speed: " + f.properties["gust_speed_ms-1"].toFixed(1) + "ms⁻¹";
    }

    // Regions
    if (props.map_style == 'regions') {
      title = "Region"
    }
    if (f.properties.NAME_1 && f.properties.NAME_0) {
      detail = f.properties.NAME_1 + ", " + f.properties.NAME_0
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
