import { Typography } from '@mui/material';
import { isNumeric, numFormat } from 'lib/helpers';

export const CompactValue = ({ label, value, maximumSignificantDigits = 3 }) => {
  if (isNumeric(value)) {
    value = numFormat(value, maximumSignificantDigits);
  }
  return (
    <div>
      <Typography variant="subtitle2" component="span">
        {label}:
      </Typography>{' '}
      <Typography variant="body2" component="span">
        {value}
      </Typography>
    </div>
  );
};
