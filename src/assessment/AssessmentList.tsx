import { PlayCircleOutline, Edit, DeleteOutline, Download } from '@mui/icons-material';
import {
  Button as IconButton,
  Paper,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  Typography,
} from '@mui/material';
import { SetterOrUpdater, useRecoilState, useSetRecoilState } from 'recoil';

import { Assessment } from 'config/assessment/assessment';
import { downloadFile } from 'lib/helpers';
import { assessmentList, currentAssessment } from 'state/assessment';
import { AssessmentUpload } from './AssessmentUpload';
import { HelpNote } from './HelpNote';
import { HeadingBox } from 'pages/HeadingBox';

function deleteAssessment(to_delete: Assessment, setAssessments: SetterOrUpdater<Assessment[]>) {
  setAssessments((prev: Assessment[]) => prev.filter((assessment) => assessment.id !== to_delete.id));
}

function downloadAssessment(assessment: Assessment) {
  const data = JSON.stringify({
    apiVersion: 1,
    assessment: assessment,
    savedAt: new Date(),
  });
  downloadFile(data, 'text/json', `assessment_${assessment.id}.json`);
}

export const AssessmentList = () => {
  const [assessments, setAssessments] = useRecoilState(assessmentList);
  const setAssessment = useSetRecoilState(currentAssessment);

  return (
    <>
      <HeadingBox>
        <Typography variant="h1">Sustainability Assessments</Typography>
      </HeadingBox>
      <div className="centred wide">
        <HelpNote sx={{ my: 2 }}>
          <p>Click "Start New" to start a sustainability assessment.</p>
          <p>Assessments are stored in the current browser window and not shared or saved online.</p>
          <p>
            To share an assessment, save it as a file using the download icon. Load shared files using the import box
            below.
          </p>
        </HelpNote>

        <TableContainer component={Paper} sx={{ my: 2 }}>
          <Table>
            <colgroup>
              <col width="10%" />
              <col width="70%" />
              <col width="20%" />
            </colgroup>
            <TableHead>
              <TableRow>
                <TableCell>#</TableCell>
                <TableCell>Assessment</TableCell>
                <TableCell>Action</TableCell>
              </TableRow>
            </TableHead>
            <TableBody>
              {assessments.map((assessment) => (
                <TableRow key={assessment.id}>
                  <TableCell>
                    <Typography
                      variant="caption"
                      style={{
                        display: 'block',
                        maxWidth: '5rem',
                        whiteSpace: 'nowrap',
                        overflow: 'hidden',
                        textOverflow: 'ellipsis',
                      }}
                    >
                      {assessment.id}
                    </Typography>
                  </TableCell>
                  <TableCell>
                    <Typography variant="body1">{assessment.description || 'Untitled'}</Typography>
                  </TableCell>
                  <TableCell>
                    <IconButton
                      onClick={() => {
                        setAssessment(assessment);
                      }}
                      title="Edit"
                      sx={{ px: 1, minWidth: '0px' }}
                    >
                      <Edit />
                    </IconButton>
                    <IconButton
                      onClick={() => {
                        deleteAssessment(assessment, setAssessments);
                      }}
                      title="Delete"
                      sx={{ px: 1, minWidth: '0px' }}
                    >
                      <DeleteOutline />
                    </IconButton>
                    <IconButton
                      onClick={() => {
                        downloadAssessment(assessment);
                      }}
                      title="Save as file"
                      sx={{ px: 1, minWidth: '0px' }}
                    >
                      <Download />
                    </IconButton>
                  </TableCell>
                </TableRow>
              ))}
              <AssessmentCreator />
            </TableBody>
          </Table>
        </TableContainer>
        <AssessmentUpload />
      </div>
    </>
  );
};

function AssessmentCreator() {
  const setAssessment = useSetRecoilState(currentAssessment);

  const addItem = () => {
    const newAssessment = new Assessment();
    setAssessment(newAssessment);
  };

  return (
    <TableRow>
      <TableCell></TableCell>
      <TableCell></TableCell>
      <TableCell>
        <IconButton variant="outlined" startIcon={<PlayCircleOutline />} onClick={addItem}>
          Start new
        </IconButton>
      </TableCell>
    </TableRow>
  );
}
