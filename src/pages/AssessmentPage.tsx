import React, { Fragment, useState } from 'react';
import _ from 'lodash';
import { TableContainer, Paper, Table, TableHead, TableRow, TableCell, TableBody, TextField, Button, Select, InputLabel, FormControl, MenuItem, Slider, Collapse, Box, IconButton } from '@mui/material';
import { Delete, KeyboardArrowDown, KeyboardArrowUp, PlayCircleOutline } from '@mui/icons-material';
import { atom, selector, useRecoilState } from 'recoil';

import { ValueLabel } from 'lib/controls/params/value-label';
import ScrollToTop from 'lib/hooks/scroll-to-top';

type IndicatorKey =
  'env_ghg'
  | 'env_air_quality'
  | 'env_energy_use'
  | 'env_habitat_disruption'
  | 'env_land'
  | 'econ_passenger'
  | 'econ_freight'
  | 'econ_passenger_occupancy'
  | 'econ_freight_load'
  | 'econ_age'
  | 'econ_road_quality'
  | 'econ_length'
  | 'econ_density'
  | 'econ_border'
  | 'soc_passenger_time'
  | 'soc_passenger_length'
  | 'soc_accidents_death'
  | 'soc_accidents_injury'
  | 'soc_accidents_death_pc'
  | 'soc_accidents_injury_pc'
  | 'soc_noise'
  | 'soc_disease'
  | 'soc_diversity'
  | 'soc_equality'
  | 'soc_inclusivity';

const INDICATOR_LABELS: ValueLabel<IndicatorKey>[] = [
  { value: 'env_ghg', label: 'GHG emissions' },
  { value: 'env_air_quality', label: 'Air quality' },
  { value: 'env_energy_use', label: 'Energy consumption (non-renewable)' },
  { value: 'env_habitat_disruption', label: 'Habitat and ecosystem disruption' },
  { value: 'env_land', label: 'Land take by transport' },

  { value: 'econ_passenger', label: 'Passenger transport volumes (passenger-km)' },
  { value: 'econ_freight', label: 'Freight transport volumes (tonne-km)' },
  { value: 'econ_passenger_occupancy', label: 'Passenger vehicle occupancy rates' },
  { value: 'econ_freight_load', label: 'Freight transport load factors' },
  { value: 'econ_age', label: 'Average age of vehicle fleet' },
  { value: 'econ_road_quality', label: 'Road quality' },
  { value: 'econ_length', label: 'Length of transport networks' },
  { value: 'econ_density', label: 'Density of transport networks' },
  { value: 'econ_border', label: 'Border restrictions' },

  { value: 'soc_passenger_time', label: 'Average passenger journey time' },
  { value: 'soc_passenger_length', label: 'Average passenger journey length' },
  { value: 'soc_accidents_death', label: 'Total number killed in traffic accidents' },
  { value: 'soc_accidents_injury', label: 'Number of injury traffic accidents' },
  { value: 'soc_accidents_death_pc', label: 'Per capita fatal accident rate' },
  { value: 'soc_accidents_injury_pc', label: 'Per capita injury accident rate' },
  { value: 'soc_noise', label: 'Population affected by traffic noise' },
  { value: 'soc_disease', label: 'Cases of respiratory disease' },
  { value: 'soc_diversity', label: 'Diversity' },
  { value: 'soc_equality', label: 'Equality and fairness' },
  { value: 'soc_inclusivity', label: 'Inclusivity' },
];

type Effect = Record<IndicatorKey, number>;
const ZERO_EFFECT = {
  'env_ghg': 0,
  'env_air_quality': 0,
  'env_energy_use': 0,
  'env_habitat_disruption': 0,
  'env_land': 0,
  'econ_passenger': 0,
  'econ_freight': 0,
  'econ_passenger_occupancy': 0,
  'econ_freight_load': 0,
  'econ_age': 0,
  'econ_road_quality': 0,
  'econ_length': 0,
  'econ_density': 0,
  'econ_border': 0,
  'soc_passenger_time': 0,
  'soc_passenger_length': 0,
  'soc_accidents_death': 0,
  'soc_accidents_injury': 0,
  'soc_accidents_death_pc': 0,
  'soc_accidents_injury_pc': 0,
  'soc_noise': 0,
  'soc_disease': 0,
  'soc_diversity': 0,
  'soc_equality': 0,
  'soc_inclusivity': 0,
}
const DEFAULT_WEIGHT = 0.5

type ScenarioKey =
  'population'
  | 'economic'
  | 'energy-cost';

type ScenarioStrength = Record<ScenarioKey, number>;

const ZERO_SCENARIO: ScenarioStrength = {
  'population': 0,
  'economic': 0,
  'energy-cost': 0
}

const SCENARIO_LABELS: ValueLabel<ScenarioKey>[] = [
  {value: 'population', label: 'Population'},
  {value: 'economic', label: 'Economic'},
  {value: 'energy-cost', label: 'Energy Costs'},
]

const SCENARIO_EFFECTS: Record<ScenarioKey, Effect> = {
  'population': {
      ...ZERO_EFFECT,
      'env_ghg': -1,
      'env_energy_use': -1,
      'econ_passenger': 1,
      'econ_freight': 1,
      'soc_accidents_death': -0.5,
      'soc_accidents_injury': -0.5,
  },
  'economic': {
      ...ZERO_EFFECT,
      'env_ghg': -1,
      'env_energy_use': -1,
      'econ_passenger': 1,
      'econ_freight': 1,
      'econ_age': 0.5,
  },
  'energy-cost': {
      ...ZERO_EFFECT,
      'env_ghg': 1,
      'env_energy_use': 1,
      'econ_passenger': -1,
      'econ_freight': -1,
  },
};

type InterventionKey =
    'infra_construction'
  | 'infra_maintenance'
  | 'demand_goods'
  | 'demand_travel'
  | 'logistics_planning'
  | 'system_eff'
  | 'fleet_eff'
  | 'fleet_elec'
  | 'road_user_charging';

type InterventionStrength = Record<InterventionKey, number>;

const ZERO_INTERVENTION: InterventionStrength = {
  'infra_construction': 0,
  'infra_maintenance': 0,
  'demand_goods': 0,
  'demand_travel': 0,
  'logistics_planning': 0,
  'system_eff': 0,
  'fleet_eff': 0,
  'fleet_elec': 0,
  'road_user_charging': 0,
}

const INTERVENTION_LABELS: ValueLabel<InterventionKey>[] = [
  {value: 'infra_construction',  label: 'Infrastructure construction'},
  {value: 'infra_maintenance',  label: 'Infrastructure maintenance'},
  {value: 'demand_goods',  label: 'Demand for goods'},
  {value: 'demand_travel',  label: 'Demand for travel'},
  {value: 'logistics_planning',  label: 'Logistics planning'},
  {value: 'system_eff',  label: 'System efficiencies'},
  {value: 'fleet_eff',  label: 'Fleet vehicle efficiencies'},
  {value: 'fleet_elec',  label: 'Fleet electrification'},
  {value: 'road_user_charging',  label: 'Road user charging'},
];

const INTERVENTION_EFFECTS: Record<InterventionKey, Effect> = {
  'infra_construction': {
      ...ZERO_EFFECT,
      'env_habitat_disruption': -1,
      'env_land': -1,
      'econ_length': 0.5,
      'econ_density': 0.5,
      'soc_passenger_length': 0.5,
  },
  'infra_maintenance': {
      ...ZERO_EFFECT,
      'env_ghg': 0.5,
      'env_energy_use': 0.5,
      'econ_road_quality': 1,
      'soc_passenger_time': 0.5,
      'soc_noise': 0.5,
  },
  'demand_goods': {
      ...ZERO_EFFECT,
      'env_ghg': -0.5,
      'env_energy_use': -0.5,
      'econ_freight': 1,
      'econ_freight_load': 1,
  },
  'demand_travel': {
      ...ZERO_EFFECT,
      'env_ghg': -0.5,
      'env_energy_use': -0.5,
      'econ_passenger': 1,
      'econ_passenger_occupancy': 1,
  },
  'logistics_planning': {
      ...ZERO_EFFECT,
      'econ_freight_load': 0.5,
      'econ_border': 0.5,
  },
  'system_eff': {
      ...ZERO_EFFECT,
      'env_ghg': 0.5,
      'env_energy_use': 0.5,
      'econ_freight': 0.5,
      'soc_passenger_time': 0.5,
  },
  'fleet_eff': {
      ...ZERO_EFFECT,
      'env_ghg': 0.5,
      'env_energy_use': 0.5,
  },
  'fleet_elec': {
    ...ZERO_EFFECT,
      'env_ghg': 1,
      'env_air_quality': 0.5,
      'soc_disease': 0.5,
  },
  'road_user_charging': {
      ...ZERO_EFFECT,
      'env_ghg': 0.5,
      'env_energy_use': 0.5,
      'econ_freight': -0.5,
  },
};

interface Assessment {
  description: string
  notes: string
  createdAt: Date

  interventionStrength: InterventionStrength
  defaultInterventionEffects: Record<InterventionKey, Effect>
  revisedInterventionEffects: Record<InterventionKey, Effect>

  scenarioStrength: ScenarioStrength
  defaultScenarioEffects: Record<ScenarioKey, Effect>
  revisedScenarioEffects: Record<ScenarioKey, Effect>

  indicatorWeights: Effect
}

const currentAssessment = atom<Assessment>({
  key: 'currentAssessment',
  default: {
    description: "",
    notes: "",
    createdAt: new Date(),

    interventionStrength: ZERO_INTERVENTION,
    defaultInterventionEffects: INTERVENTION_EFFECTS,
    revisedInterventionEffects: INTERVENTION_EFFECTS,

    scenarioStrength: ZERO_SCENARIO,
    defaultScenarioEffects: SCENARIO_EFFECTS,
    revisedScenarioEffects: SCENARIO_EFFECTS,

    indicatorWeights: _.mapValues(ZERO_EFFECT, ()=>(DEFAULT_WEIGHT)),
  }
});

const indicatorWeights = selector({
  key: 'indicatorWeights',
  get: ({get}) => {
    const assessment = get(currentAssessment);
    return assessment.indicatorWeights
  },
  set: ({get, set}, newValue) => {
    const assessment = get(currentAssessment);
    set(currentAssessment, {
      ...assessment,
      indicatorWeights: newValue
    })
  }
});

const ValueDisplay = ({value}) => {
  return (
    <Slider
      disabled
      min={-1}
      max={1}
      marks={[
        {value: -1, label: '-'},
        {value: 0, label: '0'},
        {value: 1, label: '+'},
      ]}
      track={false}
      value={value}
      />
  );
};

const WeightDisplay = ({value}) => {
  return (
    <Slider
      disabled
      min={0}
      max={1}
      marks={[
        {value: 0, label: '0'},
        {value: 1, label: '+'},
      ]}
      track={false}
      value={value}
      />
  );
};

const InterventionEffects = (
  { label, defaultEffect, revisedEffect, options, strength, setStrength, setEffect }:
    {
      label: string,
      defaultEffect: Effect,
      revisedEffect: Effect,
      options: { value: number, label: string }[],
      strength: number,
      setStrength: (e: any) => void,
      setEffect: (key: string, value: number | number[]) => void,
    }
) => {
    const [open, setOpen] = useState(false);
    return (
      <Fragment>
        <TableRow>
          <TableCell>
            <IconButton
              aria-label="expand row"
              size="small"
              onClick={() => setOpen(!open)}
            >
              {open ? <KeyboardArrowUp /> : <KeyboardArrowDown />}
            </IconButton>
          </TableCell>
          <TableCell>
            {label}
          </TableCell>
          <TableCell>
            <FormControl fullWidth sx={{ my: 1 }}>
              <InputLabel>{label}</InputLabel>
              <Select
                value={strength}
                label={label}
                onChange={setStrength}
              >
                {options.map(({ value, label }) => (
                  <MenuItem key={value} value={value}>{label}</MenuItem>
                ))}
              </Select>
              <Slider
                aria-label={label}
                value={strength}
                onChange={(e, value) => {
                  setStrength({target: {value: value}});
                }}
                step={1}
                track={false}
                marks={[
                  {value: -1, label: '-'},
                  {value: 0, label: '0'},
                  {value: 1, label: '+'},
                ]}
                min={-1}
                max={1} />
            </FormControl>
          </TableCell>
        </TableRow>
        <TableRow>
          <TableCell colSpan={3} sx={{p:0}}>
            <Collapse in={open} timeout="auto" unmountOnExit>
              <Box sx={{ m: 2 }}>
              <TableContainer component={Paper}>
                <Table size="small">
                  <TableHead>
                    <TableRow>
                      <TableCell>Indicator</TableCell>
                      <TableCell>Default Effect</TableCell>
                      <TableCell>Override</TableCell>
                      <TableCell>Current Effect</TableCell>
                    </TableRow>
                  </TableHead>
                  <TableBody>
                    {INDICATOR_LABELS.map((option) => {
                      let { value, label } = option;
                      const key = value;
                      return revisedEffect ? (
                        <TableRow key={key}>
                          <TableCell sx={{ whiteSpace: 'nowrap' }}>{label}</TableCell>
                          <TableCell>
                            <ValueDisplay value={defaultEffect[key]}/>
                          </TableCell>
                          <TableCell>
                            <Slider
                              aria-label={`${label} Revised`}
                              value={revisedEffect[key]}
                              onChange={(e, value) => {
                                setEffect(key, value);
                              } }
                              step={0.5}
                              track={false}
                              marks={[
                                {value: -1, label: '-'},
                                {value: 0, label: '0'},
                                {value: 1, label: '+'},
                              ]}
                              min={-1}
                              max={1} />
                          </TableCell>
                          <TableCell>
                            <ValueDisplay value={revisedEffect[value] * strength} />
                          </TableCell>
                        </TableRow>
                      ) : null;
                    })}
                  </TableBody>
                </Table>
              </TableContainer>
              </Box>
            </Collapse>
          </TableCell>
        </TableRow>
      </Fragment>
    );
  }

const AdjustSelector = ({ label, assessed_value, weight, setWeight}) => {
  return (
    <TableRow>
      <TableCell />
      <TableCell>
        {label}
      </TableCell>
      <TableCell>
        <ValueDisplay value={assessed_value} />
      </TableCell>
      <TableCell>
        <Slider
          aria-label={label}
          value={weight}
          onChange={setWeight}
          step={0.25}
          marks={[
            {value: 0, label: '0'},
            {value: 1, label: '+'},
          ]}
          min={0}
          max={1}
        />
      </TableCell>
      <TableCell>
      <ValueDisplay value={assessed_value * weight} />
      </TableCell>
    </TableRow>
  );
}

const IndicatorWeights = ({label, prefix, unweighted}: {
  label: string,
  prefix: string,
  unweighted: Effect,
}) => {
  const [open, setOpen] = useState(false);
  const [currentWeights, setWeights] = useRecoilState(indicatorWeights);

  // Force our way around type-checking - recoil returns the expected Effect
  // @ts-ignore
  let weights: Effect = {...currentWeights};

  let assessed_sum = 0;
  let weighted_sum = 0;
  let count = 0;
  let weight = 0;
  for (let key in unweighted) {
    if (key.includes(prefix)){
      assessed_sum += unweighted[key]
      count += 1
      weighted_sum += (unweighted[key] * weights[key])
      weight += weights[key]
    }
  }
  const assessed_value = count? assessed_sum / count: assessed_sum;
  const weighted_value = weight? weighted_sum / weight: weighted_sum;
  const total_weight = count? weight / count: count;

  return (
    <Fragment>
      <TableRow>
      <TableCell>
        <IconButton
          aria-label="expand row"
          size="small"
          onClick={() => setOpen(!open)}
          >
          {open? <KeyboardArrowUp /> : <KeyboardArrowDown />}
        </IconButton>
      </TableCell>
      <TableCell>{label}</TableCell>
      <TableCell>
        <ValueDisplay value={assessed_value} />
      </TableCell>
      <TableCell>
        <WeightDisplay value={total_weight} />
      </TableCell>
      <TableCell>
        <ValueDisplay value={weighted_value} />
      </TableCell>
    </TableRow>
    <TableRow>
      <TableCell colSpan={5} sx={{ p: 0 }}>
        <Collapse in={open} timeout="auto" unmountOnExit>
          <Table>
            <IndicatorTableColGroup/>
            <TableBody>
              {
                INDICATOR_LABELS.map((option) => {
                  let {value, label} = option;
                  const key = value;
                  if (key.includes(prefix)){
                    return (
                      <AdjustSelector
                        key={key}
                        label={label}
                        assessed_value={unweighted[key]}
                        weight={weights[key]}
                        setWeight={(_, weight) => {
                          setWeights({
                            ...weights,
                            [key]: weight
                          })
                        }}
                        />
                    )
                  }
                })
              }
            </TableBody>
          </Table>
        </Collapse>
      </TableCell>
    </TableRow>
  </Fragment>
  )
}

const IndicatorTableColGroup = () => (
  <colgroup>
    <col width="7%" />
    <col width="33%" />
    <col width="20%" />
    <col width="20%" />
    <col width="20%" />
  </colgroup>
)

export const AssessmentPage = () => {
  const [assessment, setAssessment] = useRecoilState(currentAssessment)

  let currentIndicatorsUnweighted: Effect = {...ZERO_EFFECT}
  for (let key in assessment.revisedScenarioEffects) {
    const effect: Effect = assessment.revisedScenarioEffects[key]
    const strength: number = assessment.scenarioStrength[key]
    for (let indicator in effect) {
      currentIndicatorsUnweighted[indicator] += effect[indicator] * strength
    }
  }
  for (let key in assessment.revisedInterventionEffects) {
    const effect: Effect = assessment.revisedInterventionEffects[key]
    const strength: number = assessment.interventionStrength[key]
    for (let indicator in effect) {
      currentIndicatorsUnweighted[indicator] += effect[indicator] * strength
    }
  }

  return (
    <article>
      <ScrollToTop />
{/*
      <h1>Assessments</h1>

      <TableContainer component={Paper} sx={{ my: 2 }}>
        <Table>
          <colgroup>
            <col width="10%" />
            <col width="80%" />
            <col width="10%" />
          </colgroup>
          <TableHead>
            <TableRow>
              <TableCell>#</TableCell>
              <TableCell>Assessment</TableCell>
              <TableCell>Action</TableCell>
            </TableRow>
          </TableHead>
          <TableBody>
            {
            //Map over assessments
            //- id (UUID?), name
            //- edit or delete
            //    <Button variant="outlined" startIcon={<Delete />}>
            //      Delete
            //    </Button>
            //- expand for info?
            }
            <TableRow>
              <TableCell>
                <TextField
                  disabled
                  defaultValue="1" />
              </TableCell>
              <TableCell>
                <TextField
                  required
                  fullWidth
                  label="Short description"
                  defaultValue="" />
              </TableCell>
              <TableCell>
                <Button variant="outlined" startIcon={<PlayCircleOutline />}>
                  Start
                </Button>
              </TableCell>
            </TableRow>
          </TableBody>
        </Table>
      </TableContainer>

      <br />
      <br />
      <br />
      */}

      <h2>Assessment</h2>
      <form>
        <h3>Intervention Options</h3>
        <TableContainer component={Paper} sx={{ my: 2 }}>
          <Table>
            <colgroup>
              <col width="7%" />
              <col width="43%" />
              <col width="50%" />
            </colgroup>
            <TableHead>
              <TableRow>
                <TableCell />
                <TableCell>Intervention</TableCell>
                <TableCell>Change</TableCell>
              </TableRow>
            </TableHead>
            <TableBody>
              {_.map(INTERVENTION_LABELS, (intervention) => {
                const i_key = intervention.value;
                return (
                  <InterventionEffects
                    key={intervention.value}
                    label={intervention.label}
                    defaultEffect={assessment.defaultInterventionEffects[i_key]}
                    revisedEffect={assessment.revisedInterventionEffects[i_key]}
                    strength={assessment.interventionStrength[i_key]}
                    setStrength={(e) => {
                      setAssessment({
                        ...assessment,
                        interventionStrength: {
                          ...assessment.interventionStrength,
                          [i_key]: e.target.value
                        }
                      });
                    }}
                    setEffect={(key, value) => {
                      const currentEffect = assessment.revisedInterventionEffects[i_key];

                      setAssessment({
                        ...assessment,
                        revisedInterventionEffects: {
                          ...assessment.revisedInterventionEffects,
                          [i_key]: {
                            ...currentEffect,
                            [key]: value
                          }
                        }
                      });
                    }}
                    options={[
                      { "value": -1, "label": "Decrease/Lessen" },
                      { "value": 0, "label": "No change" },
                      { "value": 1, "label": "Increase/Improve" },
                    ]} />
                );
              })}
            </TableBody>
          </Table>
        </TableContainer>

        <h3>Scenarios</h3>
        <TableContainer component={Paper} sx={{ my: 2 }}>
          <Table>
            <colgroup>
              <col width="7%" />
              <col width="43%" />
              <col width="50%" />
            </colgroup>
            <TableHead>
              <TableRow>
                <TableCell />
                <TableCell>Scenario</TableCell>
                <TableCell>Change</TableCell>
              </TableRow>
            </TableHead>
            <TableBody>
              {_.map(SCENARIO_LABELS, (intervention) => {
                const i_key = intervention.value;
                return (
                  <InterventionEffects
                    key={intervention.value}
                    label={intervention.label}
                    defaultEffect={assessment.defaultScenarioEffects[i_key]}
                    revisedEffect={assessment.revisedScenarioEffects[i_key]}
                    strength={assessment.scenarioStrength[i_key]}
                    setStrength={(e) => {
                      setAssessment({
                        ...assessment,
                        scenarioStrength: {
                          ...assessment.scenarioStrength,
                          [i_key]: e.target.value
                        }
                      });
                    }}
                    setEffect={(key, value) => {
                      const currentEffect = assessment.revisedScenarioEffects[i_key];

                      setAssessment({
                        ...assessment,
                        revisedScenarioEffects: {
                          ...assessment.revisedScenarioEffects,
                          [i_key]: {
                            ...currentEffect,
                            [key]: value
                          }
                        }
                      });
                    }}
                  options={[
                    { "value": -1, "label": "Low" },
                    { "value": 0, "label": "Central" },
                    { "value": 1, "label": "High" },
                  ]} />
                );
              })}
            </TableBody>
          </Table>
        </TableContainer>

        <h3>Impacts</h3>
        <TableContainer component={Paper} sx={{ my: 2 }}>
          <Table>
            <IndicatorTableColGroup/>
            <TableHead>
              <TableRow>
                <TableCell />
                <TableCell>Indicator</TableCell>
                <TableCell>Assessed value</TableCell>
                <TableCell>Weight</TableCell>
                <TableCell>Weighted value</TableCell>
              </TableRow>
            </TableHead>
            <TableBody>
              <IndicatorWeights
                label="Environmental"
                prefix="env_"
                unweighted={currentIndicatorsUnweighted}
                />
              <IndicatorWeights
                label="Economic"
                prefix="econ_"
                unweighted={currentIndicatorsUnweighted}
                />
              <IndicatorWeights
                label="Social"
                prefix="soc_"
                unweighted={currentIndicatorsUnweighted}
                />
            </TableBody>
          </Table>
        </TableContainer>
      </form>
      {/*
      <br />
      <br />
      <br />

      <h2>Assessment Outcomes</h2>
       */}
    </article>
  );
};
