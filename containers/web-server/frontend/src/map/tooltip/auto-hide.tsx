import { Paper } from '@mui/material';
import { createContext, useContext, useEffect, useState } from 'react';

const AutoHideContext = createContext<[hide: boolean, setHide: (empty: boolean) => void]>(null);

/**
 * Detect if the subtree has at least one desired element (marked with PreventHide)
 * Hide the container with CSS, otherwise.
 */
export const AutoHidePaper = ({ children }) => {
  const autoHideState = useState(true);
  const [hide] = autoHideState;
  return (
    <AutoHideContext.Provider value={autoHideState}>
      <Paper sx={{ display: hide ? 'none' : undefined }}>{children}</Paper>
    </AutoHideContext.Provider>
  );
};

/**
 * Use this to prevent the corresponding AutoHidePaper from being hidden
 */
export const PreventHide = () => {
  const [, setHide] = useContext(AutoHideContext);

  useEffect(() => {
    setHide(false);
  });

  return null;
};
