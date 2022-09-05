import React from 'react';
import ScrollToTop from 'lib/hooks/scroll-to-top';
import { TableContainer, Paper, Table, TableHead, TableRow, TableCell, TableBody, TextField, Button } from '@mui/material';
import Delete from '@mui/icons-material/Delete';
import PlayCircleOutline from '@mui/icons-material/PlayCircleOutline';

export const AssessmentPage = () => (
  <article>
    <ScrollToTop />
    <h1>Assessments</h1>

    <TableContainer component={Paper} sx={{ my: 2 }}>
      <Table aria-label="simple table">
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
              {/* <Button variant="outlined" startIcon={<Delete />}>
                Delete
              </Button> */}
              <Button variant="outlined" startIcon={<PlayCircleOutline />}>
                Start
              </Button>
            </TableCell>
          </TableRow>
        </TableBody>
      </Table>
    </TableContainer>
  </article>
);
