import { Paper, Table, TableBody, TableCell, TableContainer, TableFooter, TableHead, TableRow } from '@mui/material';
import { FC, ReactNode, useCallback } from 'react';
import { AssetTableRow } from './AssetTableRow';
import { ListFeature } from './use-sorted-features';

export const AssetTable: FC<{
  features: ListFeature[];
  selectedFeature?: ListFeature;
  onSelectedFeature?: (feature: ListFeature) => void;
  hoveredFeature?: ListFeature;
  onHoveredFeature?: (feature: ListFeature) => void;
  footer?: ReactNode;
}> = ({ features, selectedFeature, onSelectedFeature, hoveredFeature, onHoveredFeature, footer }) => {
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
              onMouseEnter={(feature) => onHoveredFeature?.(feature)}
              onMouseLeave={(feature) => feature === hoveredFeature && onHoveredFeature?.(null)}
            />
          ))}
        </TableBody>
        {footer && <TableFooter>{footer}</TableFooter>}
      </Table>
    </TableContainer>
  );
};
