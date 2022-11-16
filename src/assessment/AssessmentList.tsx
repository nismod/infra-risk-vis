import { PlayCircleOutline } from "@mui/icons-material";
import { Button, Paper, Table, TableBody, TableCell, TableContainer, TableHead, TableRow, TextField } from "@mui/material";

export const AssessmentList = () => (
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
  </>
);