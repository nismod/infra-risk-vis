import React, { Fragment } from 'react';
import PropTypes from 'prop-types';

const FloodControl = (props) => (
  <Fragment>
    <br />
    <h3 className="h4">Flood layers</h3>
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
        value="rcp_4p5"
        id="scenarioRadio_rcp4p5"
        onClick={(e) => props.setScenario(e.target.value)}
        />
      <label
        className="form-check-label"
        htmlFor="scenarioRadio_rcp4p5"
        >
        RCP 4.5
      </label>
    </div>
    <div className="form-check">
      <input
        className="form-check-input"
        type="radio"
        name="scenarioRadio"
        value="rcp_8p5"
        id="scenarioRadio_rcp8p5"
        onClick={(e) => props.setScenario(e.target.value)}
        />
      <label
        className="form-check-label"
        htmlFor="scenarioRadio_rcp8p5"
        >
        RCP 8.5
      </label>
    </div>

    <h4 className="h5">Flood Level</h4>
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
        &gt;4m
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
