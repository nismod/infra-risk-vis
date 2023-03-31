import { atom, useRecoilState } from 'recoil';
import { syncEffect } from 'recoil-sync';
import { date, optional } from '@recoiljs/refine';
import { Button, Dialog, DialogActions, Typography } from '@mui/material';
import { useCallback } from 'react';

const noticeAcceptedDateState = atom<Date | null>({
  key: 'noticeAcceptedDate',
  default: null,
  effects: [
    syncEffect({
      storeKey: 'local-storage',
      itemKey: 'notice-accepted',
      refine: optional(date()),
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
    <Dialog open={acceptedDate == null}>
      <Typography paragraph>
        The systemic risk analysis results shown in this tool contain licensed data that must not be shared outside the
        Government of Jamaica. By accessing the tool, you acknowledge that you understand this and agree not to download
        any data or share your access credentials with anyone else.
      </Typography>
      <DialogActions>
        <Button onClick={handleAccept}>Accept</Button>
      </DialogActions>
    </Dialog>
  );
};
