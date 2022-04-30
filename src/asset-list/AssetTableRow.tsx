import { Box, Collapse, TableCell, TableRow } from '@mui/material';
import { ListFeature } from './use-sorted-features';

export const AssetTableRow = ({
  feature,
  expanded,
  onExpandedChange,
  onMouseEnter,
  onMouseLeave,
}: {
  feature: ListFeature;
  expanded: boolean;
  onExpandedChange: (f: ListFeature, expanded: boolean) => void;
  onMouseEnter?: (f: ListFeature) => void;
  onMouseLeave?: (f: ListFeature) => void;
}) => (
  <>
    <TableRow
      onClick={() => onExpandedChange(feature, !expanded)}
      onMouseOver={() => onMouseEnter?.(feature)}
      onMouseOut={() => onMouseLeave?.(feature)}
      sx={{
        '&:hover': {
          backgroundColor: 'rgba(0, 0, 0, 0.05)',
        },
        cursor: 'pointer',
      }}
    >
      <TableCell>{feature.string_id}</TableCell>
      <TableCell>{feature.value}</TableCell>
    </TableRow>
    <TableRow>
      <TableCell style={{ paddingBottom: 0, paddingTop: 0 }} colSpan={3}>
        <Collapse in={expanded} timeout="auto" unmountOnExit>
          <Box>{feature.bbox.toString()}</Box>
        </Collapse>
      </TableCell>
    </TableRow>
  </>
);
