import { PlayCircleOutline, Edit, DeleteOutline } from "@mui/icons-material";
import { Button as IconButton, Paper, Table, TableBody, TableCell, TableContainer, TableHead, TableRow, TextField, Typography } from "@mui/material";
import { Assessment } from "config/assessment/assessment";
import { useRecoilState, useSetRecoilState } from "recoil";

import { assessmentList, currentAssessment, currentAssessmentID } from "state/assessment";

export const AssessmentList = () => {
  const [assessments, setAssessments] = useRecoilState(assessmentList);
  const setAssessmentID = useSetRecoilState(currentAssessmentID);

  function deleteAssessment(to_delete: Assessment) {
    setAssessments((prev: Assessment[]) => prev.filter((assessment) => assessment.id !== to_delete.id));
  }

  return (
    <>
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
              assessments.map((assessment) => (
                <TableRow key={assessment.id}>
                  <TableCell>
                    <Typography variant="caption" style={{display: 'block', maxWidth: '5rem', whiteSpace: 'nowrap', overflow: 'hidden', textOverflow: 'ellipsis'}}>{assessment.id}</Typography>
                  </TableCell>
                  <TableCell>
                    <Typography variant="body1">{assessment.description || "Untitled"}</Typography>
                  </TableCell>
                  <TableCell>
                    <IconButton onClick={() => {setAssessmentID(assessment.id)}} title="Edit" sx={{px: 1,minWidth: '0px'}}>
                      <Edit />
                    </IconButton>
                    <IconButton onClick={() => {deleteAssessment(assessment)}} title="Delete" sx={{px: 1,minWidth: '0px'}}>
                      <DeleteOutline />
                    </IconButton>
                  </TableCell>
                </TableRow>
              ))
            }
            <AssessmentCreator />
          </TableBody>
        </Table>
      </TableContainer>
    </>
  );
};

function AssessmentCreator() {
  const setAssessment = useSetRecoilState(currentAssessment);
  const setAssessmentID = useSetRecoilState(currentAssessmentID);

  const addItem = () => {
    const newAssessment = new Assessment();
    setAssessment(newAssessment);
    setAssessmentID(newAssessment.id);
  };

  return (
    <TableRow>
      <TableCell>
      </TableCell>
      <TableCell>
      </TableCell>
      <TableCell>
        <IconButton variant="outlined" startIcon={<PlayCircleOutline />} onClick={addItem}>
          Create
        </IconButton>
      </TableCell>
    </TableRow>
  );
}
