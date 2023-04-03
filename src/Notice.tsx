import { InfoOutlined } from '@mui/icons-material';
import { Box, Button, DialogActions, Drawer, Stack, Typography } from '@mui/material';
import { date, nullable } from '@recoiljs/refine';
import { AppLink } from 'lib/nav';
import { useCallback } from 'react';
import { atom, useRecoilState } from 'recoil';
import { syncEffect } from 'recoil-sync';

const noticeAcceptedDateState = atom<Date | null>({
  key: 'noticeAcceptedDate',
  default: null,
  effects: [
    syncEffect({
      storeKey: 'local-storage',
      itemKey: 'notice-accepted',
      refine: nullable(date()),
    }),
  ],
});

export const Notice = () => {
  const [acceptedDate, setAcceptedDate] = useRecoilState(noticeAcceptedDateState);

  const handleAccept = useCallback(
    (e) => {
      setAcceptedDate(new Date());
    },
    [setAcceptedDate],
  );

  return (
    <Drawer variant="persistent" anchor="bottom" open={acceptedDate == null}>
      <Stack direction="row" justifyContent="space-between" mx={1} my={0.5}>
        <Stack direction="row" alignItems="center" spacing={1}>
          <Box>
            <InfoOutlined color="primary" />
          </Box>
          <Typography paragraph>
            The systemic risk analysis results shown in this tool contain licensed data that must not be shared outside
            the Government of Jamaica. By accessing the tool, you acknowledge that you understand this and agree not to
            download any data or share your access credentials with anyone else.{' '}
            <AppLink to="/data">Read more about the data</AppLink>
          </Typography>
        </Stack>
        <DialogActions>
          <Button onClick={handleAccept}>Accept</Button>
        </DialogActions>
      </Stack>
    </Drawer>
  );
};
