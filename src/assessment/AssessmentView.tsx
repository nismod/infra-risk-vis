import _ from "lodash";
import { TableContainer, Paper, Table, TableHead, TableRow, TableCell, TableBody, Typography, TextField, Button } from "@mui/material";
import { useRecoilState, useSetRecoilState } from "recoil";
import { NIL as NIL_UUID } from 'uuid';

import { unweightedIndicatorSum, weightedSum } from "config/assessment/assessment";
import { Effect } from "config/assessment/effect";
import { InterventionSelection, INTERVENTION_HIERARCHY, INTERVENTION_LABELS, NO_INTERVENTIONS } from "config/assessment/interventions";
import { SCENARIO_LABELS } from "config/assessment/scenarios";
import { currentAssessment, interventionTreeConfig, interventionSelection, currentAssessmentID } from "state/assessment";

import { IndicatorTableColGroup } from "./IndicatorTableColGroup";
import { Summary } from "./Summary";
import { ValueDisplay } from "./ValueDisplay";
import { Intervention } from "./Intervention";
import { WeightGroup } from "./WeightGroup";
import { CheckboxTree } from "lib/controls/checkbox-tree/CheckboxTree";

export const AssessmentView = () => {
  const [assessment, setAssessment] = useRecoilState(currentAssessment);
  const setAssessmentID = useSetRecoilState(currentAssessmentID);
  
  const [currentInterventionsUntyped, setInterventionSelection] = useRecoilState(interventionSelection);
  // @ts-ignore: InterventionSelection
  const currentInterventions: InterventionSelection = currentInterventionsUntyped;
  let currentIndicatorsUnweighted: Effect = unweightedIndicatorSum(assessment)

  const [assessed_value, weighted_value, _total_weight] = weightedSum(
    currentIndicatorsUnweighted,
    assessment.indicatorWeights,
  );
  const createdAt = new Date(assessment.createdAt) 
  return (
    <>
      <h2>Assessment</h2>
      <Typography variant="caption" component="p">ID: {assessment.id}</Typography>
      <Typography variant="caption" component="p">Created: {createdAt.toLocaleString()}</Typography>

      <form>
        <TextField
          fullWidth
          label="Short description"
          value={assessment.description} 
          onChange={(e) => {
            setAssessment({
              ...assessment,
              description: e.target.value,
            });
          }}
          />
        <h3>Select Interventions</h3>
        <CheckboxTree
          nodes={INTERVENTION_HIERARCHY}
          config={interventionTreeConfig}
          getLabel={(node) => node.label}
          checkboxState={{checked: currentInterventions, indeterminate: NO_INTERVENTIONS}}
          onCheckboxState={(checkboxState)=>{
            // @ts-ignore: checkboxState.checked can coerce to InterventionSelection
            const nextInterventions: InterventionSelection = checkboxState.checked;
            setInterventionSelection(nextInterventions)
          }}
          expanded={[]}
          onExpanded={()=>{}}
          disableCheck={false}
          />
        <h3>Intervention Options</h3>
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
                return currentInterventions[i_key]? (
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
      <Button 
        onClick={()=>{
          setAssessmentID(NIL_UUID)
        }}>
        Save
      </Button>
    </>
  );
};