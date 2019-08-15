import React, { Fragment } from 'react';
import PropTypes from 'prop-types';

const FloodControl = (props) => (
  <Fragment>
    <h4 className="h5">Climate Scenario</h4>
    <div className="form-check">
      <input className="form-check-input" defaultChecked={true} type="radio" name="scenarioRadio" value="baseline" onClick={(e) => props.setScenario(e.target.value)}/>
      <label className="form-check-label">
        Baseline
      </label>
    </div>
    <div className="form-check">
      <input className="form-check-input" type="radio" name="scenarioRadio" value="med" onClick={(e) => props.setScenario(e.target.value)}/>
      <label className="form-check-label">
        Med
      </label>
    </div>
    <div className="form-check">
      <input className="form-check-input" type="radio" name="scenarioRadio" value="high" onClick={(e) => props.setScenario(e.target.value)}/>
      <label className="form-check-label">
        High
      </label>
    </div>

    <h4 className="h5">Flood Type</h4>
    <div className="form-check">
      <input className="form-check-input" defaultChecked={true} type="radio" name="floodtypeRadios" value="fluvial" onClick={(e) => props.setFloodType(e.target.value)}/>
      <label className="form-check-label">
        Fluvial
      </label>
    </div>
    <div className="form-check">
      <input className="form-check-input" type="radio" name="floodtypeRadios" value="pluvial" onClick={(e) => props.setFloodType(e.target.value)}/>
      <label className="form-check-label">
        Pluvial
      </label>
    </div>

    <h4 className="h5">Flood Level</h4>
    <div className="form-check">
      <input className="form-check-input" defaultChecked={true} type="checkbox" value="_50cm1m" onClick={(e) => props.setFloodLevel(e.target.value, e.target.checked)}/>
      <label className="form-check-label">
        50cm-1m
      </label>
    </div>
    <div className="form-check">
      <input className="form-check-input" defaultChecked={true} type="checkbox" value="_1m2m" onClick={(e) => props.setFloodLevel(e.target.value, e.target.checked)}/>
      <label className="form-check-label">
        1m-2m
      </label>
    </div>
    <div className="form-check">
      <input className="form-check-input" defaultChecked={true} type="checkbox" value="_2m3m" onClick={(e) => props.setFloodLevel(e.target.value, e.target.checked)}/>
      <label className="form-check-label">
        2m-3m
      </label>
    </div>
    <div className="form-check">
      <input className="form-check-input" defaultChecked={true} type="checkbox" value="_3m4m" onClick={(e) => props.setFloodLevel(e.target.value, e.target.checked)}/>
      <label className="form-check-label">
        3m-4m
      </label>
    </div>
    <div className="form-check">
      <input className="form-check-input" defaultChecked={true} type="checkbox" value="_4m999m" onClick={(e) => props.setFloodLevel(e.target.value, e.target.checked)}/>
      <label className="form-check-label">
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
