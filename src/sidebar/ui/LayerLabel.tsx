import { VisibilityOff } from '@mui/icons-material';
import { Box, Stack, Typography } from '@mui/material';
import { ShapeLegend, LegendShapeType } from 'lib/map-shapes/ShapeLegend';
import { useRecoilValue } from 'recoil';
import { zoomState } from 'state/zoom';

export const LayerLabel = ({
  label,
  type,
  color,
  minZoom,
  visible,
}: {
  label: string;
  type: LegendShapeType;
  color: string;
  minZoom?: number;
  visible: boolean;
}) => {
  return (
    <Stack direction="row" justifyContent="space-between" alignItems="center" spacing={1}>
      <Stack direction="row" alignItems="center">
        <ShapeLegend type={type} color={color} />
        <Typography>{label}</Typography>
      </Stack>
      {visible && minZoom != null && <ZoomVisibility minZoom={minZoom} />}
    </Stack>
  );
};

function ZoomVisibility({ minZoom }: { minZoom: number }) {
  const currentZoom = useRecoilValue(zoomState);

  return (
    currentZoom < minZoom && (
      <Box title="This data is only visible at higher map zoom levels">
        <VisibilityOff fontSize="small" color="disabled" />
      </Box>
    )
  );
}
