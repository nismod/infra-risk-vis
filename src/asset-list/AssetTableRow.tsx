import { Box, Collapse, TableCell, TableRow } from '@mui/material';
import { ListFeature } from 'pages/AssetListPage';

export const AssetTableRow = ({
  feature,
  expanded,
  onExpandedChange,
}: {
  feature: ListFeature;
  expanded: boolean;
  onExpandedChange: (f: ListFeature, expanded: boolean) => void;
}) => (
  <>
    <TableRow
      onClick={() => onExpandedChange(feature, !expanded)}
      sx={{
        '&:hover': {
          backgroundColor: 'rgba(0, 0, 0, 0.05)',
        },
        cursor: 'pointer',
      }}
    >
      <TableCell>{feature.id}</TableCell>
      <TableCell>{feature.string_id}</TableCell>
      <TableCell>{feature.value}</TableCell>
    </TableRow>
    <TableRow>
      <TableCell style={{ paddingBottom: 0, paddingTop: 0 }} colSpan={3}>
        <Collapse in={expanded} timeout="auto" unmountOnExit>
          <Box>{feature.bbox_wkt}</Box>
        </Collapse>
      </TableCell>
    </TableRow>
  </>
);
