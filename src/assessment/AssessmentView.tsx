import _ from 'lodash';
import {
  TableContainer,
  Paper,
  Table,
  TableHead,
  TableRow,
  TableCell,
  TableBody,
  Typography,
  TextField,
  Button,
  Stack,
  FormControlLabel,
  Checkbox,
} from '@mui/material';
import { useRecoilState, useSetRecoilState } from 'recoil';

import { unweightedIndicatorSum, weightedSum } from 'config/assessment/assessment';
import { Effect } from 'config/assessment/effect';
import { InterventionSelection, INTERVENTION_LABELS } from 'config/assessment/interventions';
import { SCENARIO_LABELS } from 'config/assessment/scenarios';
import { currentAssessment, interventionSelection, currentAssessmentInList } from 'state/assessment';

import { IndicatorTableColGroup } from './IndicatorTableColGroup';
import { Summary } from './Summary';
import { ValueDisplay } from './ValueDisplay';
import { Intervention } from './Intervention';
import { WeightGroup } from './WeightGroup';
import { HelpNote } from './HelpNote';

export const AssessmentView = () => {
  const [assessment, setAssessment] = useRecoilState(currentAssessment);
  const setAssessmentInList = useSetRecoilState(currentAssessmentInList);

  const [currentInterventionsUntyped, setInterventionSelection] = useRecoilState(interventionSelection);
  // @ts-ignore: InterventionSelection
  const currentInterventions: InterventionSelection = currentInterventionsUntyped;
  let currentIndicatorsUnweighted: Effect = unweightedIndicatorSum(assessment);

  const [assessed_value, weighted_value, _total_weight] = weightedSum(
    currentIndicatorsUnweighted,
    assessment.indicatorWeights,
  );
  const createdAt = new Date(assessment.createdAt);
  return (
    <>
      <Typography variant="caption" component="p">
        ID: {assessment.id}
      </Typography>
      <Typography variant="caption" component="p" sx={{ mb: 2 }}>
        Created: {createdAt.toLocaleString()}
      </Typography>
      <h2>Sustainability Assessment</h2>

      <form>
        <TextField
          fullWidth
          sx={{ my: 1 }}
          label="Title"
          value={assessment.description}
          onChange={(e) => {
            setAssessment({
              ...assessment,
              description: e.target.value,
            });
          }}
        />
        <TextField
          fullWidth
          sx={{ my: 1 }}
          multiline
          minRows={3}
          label="Notes"
          value={assessment.notes}
          onChange={(e) => {
            setAssessment({
              ...assessment,
              notes: e.target.value,
            });
          }}
        />
        <h3>Select Interventions</h3>
        <HelpNote>
          <p>
            Select one or more interventions for assessment. This will provide a template with some preset values for
            the sustainability indicators.
          </p>
          <p>If none of the below are relevant, select "Custom intervention".</p>
        </HelpNote>
        <TableContainer component={Paper} sx={{ my: 2, px: 1 }}>
          <Table>
            <colgroup>
              <col width="30%" />
              <col width="70%" />
            </colgroup>
            <TableHead>
              <TableRow>
                <TableCell>Intervention</TableCell>
                <TableCell>Description</TableCell>
              </TableRow>
            </TableHead>
            <TableBody>
              {INTERVENTION_LABELS.map((intervention) => (
                <TableRow>
                  <TableCell>
                    <FormControlLabel
                      key={intervention.value}
                      label={intervention.label}
                      style={{ width: '100%' }}
                      control={
                        <Checkbox
                          checked={currentInterventions[intervention.value]}
                          onChange={(event) => {
                            const nextInterventions: InterventionSelection = { ...currentInterventions };
                            nextInterventions[intervention.value] = event.currentTarget.checked;
                            setInterventionSelection(nextInterventions);
                          }}
                          onClick={(e) => e.stopPropagation()}
                          style={{ pointerEvents: 'auto' }}
                        />
                      }
                    ></FormControlLabel>
                  </TableCell>
                  <TableCell>{intervention.description}</TableCell>
                </TableRow>
              ))}
            </TableBody>
          </Table>
        </TableContainer>
        <h3>Intervention Options</h3>
        <HelpNote>
          <p>The table below shows the interventions that have been selected for evaluation.</p>
          <p>
            Choose a direction of change for each intervention. The default is "No intervention" and assumes doing
            nothing. For example, set this to "Increase/Improve" in order to evaluate a project or positive change.
          </p>
          <p>
            Expand each row to evaluate the impact of the intervention against each indicator. This will show the
            default (preset) value and any changes to the indicators. Make a note of reasons for any changes from the
            default values.
          </p>
        </HelpNote>
        <TableContainer component={Paper} sx={{ my: 2, px: 1 }}>
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
                return currentInterventions[i_key] ? (
                  <Intervention
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
                          [i_key]: e.target.value,
                        },
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
                            [key]: value,
                          },
                        },
                      });
                    }}
                    options={[
                      { value: -1, label: 'Decrease/Lessen' },
                      { value: 0, label: 'No intervention' },
                      { value: 1, label: 'Increase/Improve' },
                    ]}
                  />
                ) : null;
              })}
            </TableBody>
          </Table>
        </TableContainer>

        <h3>Scenarios</h3>
        <HelpNote>
          <p>The table below shows the scenarios that have been selected for evaluation.</p>
          <p>
            Choose a direction of change for each scenario. The default is "Central" which gives a baseline scenario
            with neutral effects. For example, test interventions under a range of scenarios to evaluate how there may
            be a range of sustainability outcomes.
          </p>
          <p>
            Expand each row to evaluate the impact of the scenario against each indicator. This will show the default
            (preset) value and any changes to the indicators. Make a note of reasons for any changes from the default
            values.
          </p>
        </HelpNote>
        <TableContainer component={Paper} sx={{ my: 2, px: 1 }}>
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
                  <Intervention
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
                          [i_key]: e.target.value,
                        },
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
                            [key]: value,
                          },
                        },
                      });
                    }}
                    options={[
                      { value: -1, label: 'Low' },
                      { value: 0, label: 'Central' },
                      { value: 1, label: 'High' },
                    ]}
                  />
                );
              })}
            </TableBody>
          </Table>
        </TableContainer>

        <h3>Impacts</h3>
        <HelpNote>
          <p>
            The table below gives a summary of the indicator values chosen, averaging the effects of interventions and
            scenarios, then grouping into environmental, economic, and social sustainability.
          </p>
          <p>
            Expand each row to change the weighting given to each indicator. By default, each indicator is weighted
            equally at 0.5.
          </p>
        </HelpNote>
        <TableContainer component={Paper} sx={{ my: 2 }}>
          <Table>
            <IndicatorTableColGroup />
            <TableHead>
              <TableRow className="group-header">
                <TableCell />
                <TableCell>Indicator</TableCell>
                <TableCell>Assessed value</TableCell>
                <TableCell>Weight</TableCell>
                <TableCell>Weighted value</TableCell>
              </TableRow>
            </TableHead>
            <TableBody>
              <WeightGroup label="Environmental" prefix="env" unweighted={currentIndicatorsUnweighted} />
              <WeightGroup label="Economic" prefix="econ" unweighted={currentIndicatorsUnweighted} />
              <WeightGroup label="Social" prefix="soc" unweighted={currentIndicatorsUnweighted} />
              <TableRow className="group-header">
                <TableCell />
                <TableCell component="th">
                  <strong style={{ fontWeight: 500 }}>Overall</strong>
                </TableCell>
                <TableCell>
                  <ValueDisplay value={assessed_value} />
                </TableCell>
                <TableCell></TableCell>
                <TableCell>
                  <ValueDisplay value={weighted_value} />
                </TableCell>
              </TableRow>
            </TableBody>
          </Table>
        </TableContainer>
      </form>
      <Summary
        overall_assessed={assessed_value}
        overall_weighted={weighted_value}
        unweighted={currentIndicatorsUnweighted}
      />
      <Stack direction="row" justifyContent="flex-end" alignItems="flex-start" spacing={2}>
        <Button
          variant="contained"
          color="primary"
          onClick={() => {
            setAssessmentInList(assessment);
            setAssessment(undefined);
          }}
        >
          Save
        </Button>
        <Button
          variant="outlined"
          color="error"
          onClick={() => {
            setAssessment(undefined);
          }}
        >
          Discard
        </Button>
      </Stack>
    </>
  );
};
