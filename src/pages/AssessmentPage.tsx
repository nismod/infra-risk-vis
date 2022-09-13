import React, { Fragment, useState } from 'react';
import _ from 'lodash';
import { TableContainer, Paper, Table, TableHead, TableRow, TableCell, TableBody, TextField, Button, Select, InputLabel, FormControl, MenuItem, Slider, Collapse, Box, IconButton } from '@mui/material';
import { Delete, KeyboardArrowDown, KeyboardArrowUp, PlayCircleOutline } from '@mui/icons-material';
import { atom, useRecoilState, useRecoilValue } from 'recoil';

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

type ScenarioKey =
  'population'
  | 'economic'
  | 'energy-cost';

type InterventionKey =
  'fleet_elec'
  | 'fleet_eff'
  | 'system_eff'
  | 'demand_goods'
  | 'demand_travel'
  | 'infra_construction'
  | 'infra_maintenance'
  | 'logistics_planning'
  | 'road_user_charging';

class Intervention {
  key: string
  label: string
  effect: Effect
}

const SCENARIOS: Record<ScenarioKey, Intervention> = {
  'population': {
    key: 'population',
    label: 'Population',
    effect: {
      ...ZERO_EFFECT,
      'env_ghg': -1,
      'env_energy_use': -1,
      'econ_passenger': 1,
      'econ_freight': 1,
      'soc_accidents_death': -0.5,
      'soc_accidents_injury': -0.5,
    }
  },
  'economic': {
    key: 'economic',
    label: 'Economic',
    effect: {
      ...ZERO_EFFECT,
      'env_ghg': -1,
      'env_energy_use': -1,
      'econ_passenger': 1,
      'econ_freight': 1,
      'econ_age': 0.5,
    }
  },
  'energy-cost': {
    key: 'energy-cost',
    label: 'Energy',
    effect: {
      ...ZERO_EFFECT,
      'env_ghg': 1,
      'env_energy_use': 1,
      'econ_passenger': -1,
      'econ_freight': -1,
    }
  },
};

const INTERVENTIONS: Record<InterventionKey, Intervention> = {
  'fleet_elec': {
    key: 'fleet_elec',
    label: 'Fleet electrification',
    effect: {
      ...ZERO_EFFECT,
      'env_ghg': 1,
      'env_air_quality': 0.5,
      'soc_disease': 0.5,
    }
  },
  'fleet_eff': {
    key: 'fleet_eff',
    label: 'Fleet vehicle efficiencies',
    effect: {
      ...ZERO_EFFECT,
      'env_ghg': 0.5,
      'env_energy_use': 0.5,
    }
  },
  'system_eff': {
    key: 'system_eff',
    label: 'System efficiencies',
    effect: {
      ...ZERO_EFFECT,
      'env_ghg': 0.5,
      'env_energy_use': 0.5,
      'econ_freight': 0.5,
      'soc_passenger_time': 0.5,
    }
  },
  'demand_goods': {
    key: 'demand_goods',
    label: 'Demand for goods',
    effect: {
      ...ZERO_EFFECT,
      'env_ghg': -0.5,
      'env_energy_use': -0.5,
      'econ_freight': 1,
      'econ_freight_load': 1,
    }
  },
  'demand_travel': {
    key: 'demand_travel',
    label: 'Demand for travel',
    effect: {
      ...ZERO_EFFECT,
      'env_ghg': -0.5,
      'env_energy_use': -0.5,
      'econ_passenger': 1,
      'econ_passenger_occupancy': 1,
    }
  },
  'infra_construction': {
    key: 'infra_construction',
    label: 'Infrastructure construction',
    effect: {
      ...ZERO_EFFECT,
      'env_habitat_disruption': -1,
      'env_land': -1,
      'econ_length': 0.5,
      'econ_density': 0.5,
      'soc_passenger_length': 0.5,
    }
  },
  'infra_maintenance': {
    key: 'infra_maintenance',
    label: 'Infrastructure maintenance',
    effect: {
      ...ZERO_EFFECT,
      'env_ghg': 0.5,
      'env_energy_use': 0.5,
      'econ_road_quality': 1,
      'soc_passenger_time': 0.5,
      'soc_noise': 0.5,
    }
  },
  'logistics_planning': {
    key: 'logistics_planning',
    label: 'Logistics planning',
    effect: {
      ...ZERO_EFFECT,
      'econ_freight_load': 0.5,
      'econ_border': 0.5,
    }
  },
  'road_user_charging': {
    key: 'road_user_charging',
    label: 'Road user charging',
    effect: {
      ...ZERO_EFFECT,
      'env_ghg': 0.5,
      'env_energy_use': 0.5,
      'econ_freight': -0.5,
    }
  },
};

class Assessment {
  description: string = ""
  notes: string = ""
  createdAt: Date = new Date()
  scenarios: Record<ScenarioKey, Intervention> = {...SCENARIOS}
  interventions: Record<InterventionKey, Intervention> = {...INTERVENTIONS}
  getScenarios() {
    let scenarios = []
    for (let key in this.scenarios){
      scenarios.push(this.scenarios[key])
    }
    return scenarios
  }
  getInterventions() {
    let interventions = []
    for (let key in this.interventions){
      interventions.push(this.interventions[key])
    }
    return interventions
  }
}

const currentAssessment = atom<Assessment>({
  key: 'currentAssessment',
  default: new Assessment()
});

type InterventionStrength = Record<InterventionKey, number>;
const interventionStrength = atom<InterventionStrength>({
  key: 'interventionStrength',
  default: {
    'fleet_elec': 0,
    'fleet_eff': 0,
    'system_eff': 0,
    'demand_goods': 0,
    'demand_travel': 0,
    'infra_construction': 0,
    'infra_maintenance': 0,
    'logistics_planning': 0,
    'road_user_charging': 0,
  }
});
type ScenarioStrength = Record<ScenarioKey, number>;
const scenarioStrength = atom<ScenarioStrength>({
  key: 'scenarioStrength',
  default: {
    'population': 0,
    'economic': 0,
    'energy-cost': 0
  }
});
type IndicatorWeight = Record<IndicatorKey, number>;
const indicatorWeights = atom<IndicatorWeight>({
  key: 'indicatorWeights',
  default: {
    'env_ghg': 1,
    'env_air_quality': 1,
    'env_energy_use': 1,
    'env_habitat_disruption': 1,
    'env_land': 1,
    'econ_passenger': 1,
    'econ_freight': 1,
    'econ_passenger_occupancy': 1,
    'econ_freight_load': 1,
    'econ_age': 1,
    'econ_road_quality': 1,
    'econ_length': 1,
    'econ_density': 1,
    'econ_border': 1,
    'soc_passenger_time': 1,
    'soc_passenger_length': 1,
    'soc_accidents_death': 1,
    'soc_accidents_injury': 1,
    'soc_accidents_death_pc': 1,
    'soc_accidents_injury_pc': 1,
    'soc_noise': 1,
    'soc_disease': 1,
    'soc_diversity': 1,
    'soc_equality': 1,
    'soc_inclusivity': 1,
  }
})
const revisedInterventions = atom<Record<InterventionKey, Intervention>>({
  key: 'revisedInterventions',
  default: INTERVENTIONS
})
const revisedScenarios = atom<Record<ScenarioKey, Intervention>>({
  key: 'revisedScenarios',
  default: SCENARIOS
})

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
                          <TableCell>{defaultEffect[key]}</TableCell>
                          <TableCell>
                            <Slider
                              aria-label={`${label} Revised`}
                              value={revisedEffect[key]}
                              onChange={(e, value) => {
                                setEffect(key, value);
                              } }
                              step={0.5}
                              track={false}
                              marks
                              min={-1}
                              max={1} />
                          </TableCell>
                          <TableCell>{revisedEffect[value] * strength}</TableCell>
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
        {assessed_value}
      </TableCell>
      <TableCell>
        <Slider
          aria-label={label}
          value={weight}
          onChange={setWeight}
          step={0.5}
          marks
          min={0}
          max={2}
        />
      </TableCell>
      <TableCell>
        {assessed_value * weight}
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

  let assessed_value = 0;
  let weighted_value = 0;
  for (let key in unweighted) {
    if (key.includes(prefix)){
      assessed_value += unweighted[key]
      weighted_value += (unweighted[key] * currentWeights[key])
    }
  }

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
      <TableCell>{assessed_value}</TableCell>
      <TableCell>-</TableCell>
      <TableCell>{weighted_value}</TableCell>
    </TableRow>
    <TableRow>
      <TableCell colSpan={5} sx={{ p: 0 }}>
        <Collapse in={open} timeout="auto" unmountOnExit>
          <Box sx={{ pl: 0 }}>
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
                          weight={currentWeights[key]}
                          setWeight={(_, weight) => {
                            setWeights({...currentWeights, [key]: weight})
                          }}
                          />
                      )
                    }
                  })
                }
              </TableBody>
            </Table>
          </Box>
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
  let assessment = useRecoilValue(currentAssessment)
  const [currentInterventionStrength, setInterventionStrength] = useRecoilState(interventionStrength);
  const [currentScenarioStrength, setScenarioStrength] = useRecoilState(scenarioStrength);
  const [currentRevisedScenarios, setRevisedScenarios] = useRecoilState(revisedScenarios);
  const [currentRevisedInterventions, setRevisedInterventions] = useRecoilState(revisedInterventions);

  let currentIndicatorsUnweighted: Effect = {...ZERO_EFFECT}
  for (let key in currentRevisedScenarios) {
    const effect: Effect = currentRevisedScenarios[key].effect
    const strength: number = currentScenarioStrength[key]
    for (let indicator in effect) {
      currentIndicatorsUnweighted[indicator] += effect[indicator] * strength
    }
  }
  for (let key in currentRevisedInterventions) {
    const effect: Effect = currentRevisedInterventions[key].effect
    const strength: number = currentInterventionStrength[key]
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
              {assessment.getInterventions().map((intervention) => (
                <InterventionEffects
                  key={intervention.key}
                  label={intervention.label}
                  defaultEffect={intervention.effect}
                  revisedEffect={currentRevisedInterventions[intervention.key].effect}
                  strength={currentInterventionStrength[intervention.key]}
                  setStrength={(e) => {
                    setInterventionStrength({...currentInterventionStrength, [intervention.key]: e.target.value})
                  }}
                  setEffect={(key, value) => {
                    const currentEffect = currentRevisedInterventions[intervention.key].effect;

                    setRevisedInterventions({
                      ...currentRevisedInterventions,
                      [intervention.key]: {
                        ...intervention,
                        effect: {
                          ...currentEffect,
                          [key]: value
                        }
                      }
                    })
                  }}
                  options={[
                    { "value": -1, "label": "Decrease/Lessen" },
                    { "value": 0, "label": "No change" },
                    { "value": 1, "label": "Increase/Improve" },
                  ]} />
              ))}
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
              {assessment.getScenarios().map((intervention) => (
                <InterventionEffects
                  key={intervention.key}
                  label={intervention.label}
                  defaultEffect={intervention.effect}
                  revisedEffect={currentRevisedScenarios[intervention.key].effect}
                  strength={currentScenarioStrength[intervention.key]}
                  setStrength={(e) => {
                    setScenarioStrength({...currentScenarioStrength, [intervention.key]: e.target.value})
                  }}
                  setEffect={(key, value) => {
                    const currentEffect = currentRevisedScenarios[intervention.key].effect;

                    setRevisedScenarios({
                      ...currentRevisedScenarios,
                      [intervention.key]: {
                        ...intervention,
                        effect: {
                          ...currentEffect,
                          [key]: value
                        }
                      }
                    })
                  }}
                  options={[
                    { "value": -1, "label": "Low" },
                    { "value": 0, "label": "Central" },
                    { "value": 1, "label": "High" },
                  ]} />
              ))}
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
