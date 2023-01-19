import { Slider as BaseSlider } from '@mui/material';
import { styled } from '@mui/material/styles';

export const Slider = styled(BaseSlider)(({ theme }) => ({
  marginTop: '15px',
  '& .MuiSlider-valueLabel': {
    fontSize: 12,
    fontWeight: 'normal',
    top: 0,
    backgroundColor: 'unset',
    color: theme.palette.text.primary,
    '&:before': {
      display: 'none',
    },
    '& *': {
      background: 'transparent',
      color: theme.palette.mode === 'dark' ? '#fff' : '#000',
    },
  },
}));

export const ValueDisplay = ({ value }: { value: number }) => (
  <>
    <Slider
      disabled
      min={-1}
      max={1}
      marks={[
        { value: -1, label: '-' },
        { value: 0, label: '0' },
        { value: 1, label: '+' },
      ]}
      track={false}
      value={value}
      valueLabelDisplay="on"
    />
  </>
);
