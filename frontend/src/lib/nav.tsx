import { LinkProps, Link as MuiLink } from '@mui/material';
import { Link as RouterLink } from 'react-router-dom';

export const ExtLink = ({ ...props }: Omit<LinkProps<'a'>, 'component'>) => {
  return <MuiLink target="_blank" rel="noopener noreferrer" {...props} />;
};

export const AppLink = ({ ...props }: Omit<LinkProps<RouterLink>, 'component'>) => {
  return <MuiLink component={RouterLink} {...props} />;
};
