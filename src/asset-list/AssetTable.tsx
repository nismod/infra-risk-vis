import { Paper, Table, TableBody, TableCell, TableContainer, TableFooter, TableHead, TableRow } from '@mui/material';
import { ListFeature } from 'pages/AssetListPage';
import { useCallback } from 'react';
import { AssetTableRow } from './AssetTableRow';

export const AssetTable = ({ features, selectedFeature, onSelectedFeature, footer }) => {
  const handleExpandedChange = useCallback(
    (feature: ListFeature, expanded: boolean) => {
      onSelectedFeature(expanded ? feature : null);
    },
    [onSelectedFeature],
  );

  return (
    <TableContainer component={Paper}>
      <Table stickyHeader size="small">
        <TableHead>
          <TableRow>
            <TableCell width={100}>ID</TableCell>
            <TableCell width={150}>Text ID</TableCell>
            <TableCell>Value</TableCell>
          </TableRow>
        </TableHead>
        <TableBody sx={{ width: '100%', maxHeight: '50vh', overflowY: 'scroll' }}>
          {features.map((feature) => (
            <AssetTableRow
              key={feature.id}
              feature={feature}
              expanded={feature.id === selectedFeature?.id}
              onExpandedChange={handleExpandedChange}
            />
          ))}
        </TableBody>
        {footer && <TableFooter>{footer}</TableFooter>}
      </Table>
    </TableContainer>
  );
};
