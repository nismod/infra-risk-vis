import React, { Fragment } from 'react';
import PropTypes from 'prop-types';

const FloodControl = (props) => (
  <Fragment>
    <h4 className="h5">Climate Scenario</h4>
    <div className="form-check">
      <input
        className="form-check-input"
        defaultChecked={true} type="radio"
        name="scenarioRadio"
        value="baseline"
        id="scenarioRadio_baseline"
        onClick={(e) => props.setScenario(e.target.value)}
        />
      <label
        className="form-check-label"
        htmlFor="scenarioRadio_baseline"
        >
        Baseline
      </label>
    </div>
    <div className="form-check">
      <input
        className="form-check-input"
        type="radio"
        name="scenarioRadio"
        value="med"
        id="scenarioRadio_med"
        onClick={(e) => props.setScenario(e.target.value)}
        />
      <label
        className="form-check-label"
        htmlFor="scenarioRadio_med"
        >
        Med
      </label>
    </div>
    <div className="form-check">
      <input
        className="form-check-input"
        type="radio"
        name="scenarioRadio"
        value="high"
        id="scenarioRadio_high"
        onClick={(e) => props.setScenario(e.target.value)}
        />
      <label
        className="form-check-label"
        htmlFor="scenarioRadio_high"
        >
        High
      </label>
    </div>

    <h4 className="h5">Flood Type</h4>
    <div className="form-check">
      <input
        className="form-check-input"
        defaultChecked={true} type="radio"
        name="floodtypeRadios"
        value="fluvial"
        id="floodtypeRadios_fluvial"
        onClick={(e) => props.setFloodType(e.target.value)}
        />
      <label
        className="form-check-label"
        htmlFor="floodtypeRadios_fluvial"
        >
        Fluvial
      </label>
    </div>
    <div className="form-check">
      <input
        className="form-check-input"
        type="radio"
        name="floodtypeRadios"
        value="pluvial"
        id="floodtypeRadios_pluvial"
        onClick={(e) => props.setFloodType(e.target.value)}
        />
      <label
        className="form-check-label"
        htmlFor="floodtypeRadios_pluvial"
        >
        Pluvial
      </label>
    </div>

    <h4 className="h5">Flood Level</h4>
    <div className="form-check">
      <input
        className="form-check-input"
        defaultChecked={true}
        type="checkbox"
        name="floodLevelCheck"
        value="_50cm1m"
        id="floodLevelCheck_50cm1m"
        onClick={(e) => props.setFloodLevel(e.target.value, e.target.checked)}
        />
      <span
        className="dot"
        style={{backgroundColor: "#ffffff", boxShadow: "black 0px 0px 1px 1px"}}>
      </span>
      <label
        className="form-check-label"
        htmlFor="floodLevelCheck_50cm1m"
        >
        50cm-1m
      </label>
    </div>
    <div className="form-check">
      <input
        className="form-check-input"
        defaultChecked={true}
        type="checkbox"
        name="floodLevelCheck"
        value="_1m2m"
        id="floodLevelCheck_1m2m"
        onClick={(e) => props.setFloodLevel(e.target.value, e.target.checked)}
        />
      <span
        className="dot"
        style={{backgroundColor: "#58cced"}}>
      </span>
      <label
        className="form-check-label"
        htmlFor="floodLevelCheck_1m2m"
        >
        1m-2m
      </label>
    </div>
    <div className="form-check">
      <input
        className="form-check-input"
        defaultChecked={true}
        type="checkbox"
        name="floodLevelCheck"
        value="_2m3m"
        id="floodLevelCheck_2m3m"
        onClick={(e) => props.setFloodLevel(e.target.value, e.target.checked)}
        />
      <span
        className="dot"
        style={{backgroundColor: "#3895d3"}}>
      </span>
      <label
        className="form-check-label"
        htmlFor="floodLevelCheck_2m3m"
        >
        2m-3m
      </label>
    </div>
    <div className="form-check">
      <input
        className="form-check-input"
        defaultChecked={true}
        type="checkbox"
        name="floodLevelCheck"
        value="_3m4m"
        id="floodLevelCheck_3m4m"
        onClick={(e) => props.setFloodLevel(e.target.value, e.target.checked)}
        />
      <span
        className="dot"
        style={{backgroundColor: "#1261a0"}}>
      </span>
      <label
        className="form-check-label"
        htmlFor="floodLevelCheck_3m4m"
        >
        3m-4m
      </label>
    </div>
    <div className="form-check">
      <input
        className="form-check-input"
        defaultChecked={true}
        type="checkbox"
        name="floodLevelCheck"
        value="_4m999m"
        id="floodLevelCheck_4m999m"
        onClick={(e) => props.setFloodLevel(e.target.value, e.target.checked)}
        />
      <span
        className="dot"
        style={{backgroundColor: "#072f5f"}}>
      </span>
      <label
        className="form-check-label"
        htmlFor="floodLevelCheck_4m999m"
        >
        >4m
      </label>
    </div>
  </Fragment>
)

FloodControl.propTypes = {
  setScenario: PropTypes.func,
  setFloodType: PropTypes.func,
  setFloodLevel: PropTypes.func
}

export default FloodControl;
