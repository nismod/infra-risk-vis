import { List } from '@mui/material';

import { DetailHeader, IdSubheader } from '@/details/features/detail-components';

export const WdpaDetails = ({ f }) => {
  return (
    <>
      <DetailHeader>{f.NAME}</DetailHeader>
      <IdSubheader id={f.WDPA_PID} />
      <List>{/* TODO */}</List>
    </>
  );
};
