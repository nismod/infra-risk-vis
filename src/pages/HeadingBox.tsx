import { Paper } from '@mui/material';
import { styled } from '@mui/material/styles';

export const HeadingBox = styled(Paper)(({ theme }) => ({
  backgroundColor: theme.palette.primary.main,
  color: '#fff',
  paddingTop: '8rem',

  [theme.breakpoints.up('sm')]: {
    paddingTop: '16rem',
  },
  paddingLeft: theme.spacing(4),
  paddingRight: theme.spacing(4),
  paddingBottom: theme.spacing(2),
  borderRadius: 0,
}));
