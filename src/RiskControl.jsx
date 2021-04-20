import React, { Fragment } from 'react';

const labels = {
  EAEL: "Expected Annual Economic Losses (EAEL, million US$)",
  EAD: "Max Expected Annual Damages (EAD, million US$)",
  total: "Max Total Expected Risk (Max EAD + EAEL for 30 day disruption, million US$)"
}

const RiskControl = (props) => (
  <Fragment>
  <br />
  <h3 className="h4">Select risk metric</h3>
  <div className="form-check">
    <input
      className="form-check-input"
      defaultChecked={props.riskMetric === 'total'} type="radio"
      name="riskRadio"
      value="total"
      id="riskRadio_total"
      onClick={(e) => props.setRiskMetric(e.target.value)}
      />
    <label
      className="form-check-label"
      htmlFor="riskRadio_total"
      >
      {labels.total}
    </label>
  </div>
  <div className="form-check">
    <input
      className="form-check-input"
      defaultChecked={props.riskMetric === 'EAD'} type="radio"
      name="riskRadio"
      value="EAD"
      id="riskRadio_EAD"
      onClick={(e) => props.setRiskMetric(e.target.value)}
      />
    <label
      className="form-check-label"
      htmlFor="riskRadio_EAD"
      >
      {labels.EAD}
    </label>
  </div>
  <div className="form-check">
    <input
      className="form-check-input"
      defaultChecked={props.riskMetric === 'EAEL'} type="radio"
      name="riskRadio"
      value="EAEL"
      id="riskRadio_EAEL"
      onClick={(e) => props.setRiskMetric(e.target.value)}
      />
    <label
      className="form-check-label"
      htmlFor="riskRadio_EAEL"
      >
      {labels.EAEL}
    </label>
  </div>

  <small>{labels[props.riskMetric]}</small>
  <svg width="270" height="25" version="1.1" xmlns="http://www.w3.org/2000/svg">
    <defs>
      <linearGradient id="gradient" x1="0" x2="1" y1="0" y2="0">
        <stop offset="0%" stopColor="#fcfcb8" />
        <stop offset="20%" stopColor="#ff9c66" />
        <stop offset="40%" stopColor="#d03f6f" />
        <stop offset="60%" stopColor="#792283" />
        <stop offset="80%" stopColor="#3f0a72" />
        <stop offset="100%" stopColor="#151030" />
      </linearGradient>
    </defs>
    <g fill="none" fontSize="10" fontFamily="sans-serif">
    </g>
    <rect x="2" y="0" width="258" height="10" fill="url(#gradient)"/>
    <g fill="none" fontSize="10" transform="translate(2,10)" fontFamily="sans-serif" textAnchor="middle">
      <g transform="translate(0.5,0)">
        <line stroke="currentColor" y2="3"></line>
        <text fill="currentColor" y="6" dy="0.71em">0</text>
      </g>
      <g transform="translate(43,0)">
        <line stroke="currentColor" y2="3"></line>
        <text fill="currentColor" y="6" dy="0.71em">0.001</text>
      </g>
      <g transform="translate(86,0)">
        <line stroke="currentColor" y2="3"></line>
        <text fill="currentColor" y="6" dy="0.71em">0.01</text>
      </g>
      <g transform="translate(129,0)">
        <line stroke="currentColor" y2="3"></line>
        <text fill="currentColor" y="6" dy="0.71em">0.1</text>
      </g>
      <g transform="translate(172,0)">
        <line stroke="currentColor" y2="3"></line>
        <text fill="currentColor" y="6" dy="0.71em">1</text>
      </g>
      <g transform="translate(215,0)">
        <line stroke="currentColor" y2="3"></line>
        <text fill="currentColor" y="6" dy="0.71em">10</text>
      </g>
      <g transform="translate(257.5,0)">
        <line stroke="currentColor" y2="3"></line>
        <text fill="currentColor" y="6" dy="0.71em">100</text>
      </g>
    </g>
  </svg>
  </Fragment>
);

export default RiskControl;
