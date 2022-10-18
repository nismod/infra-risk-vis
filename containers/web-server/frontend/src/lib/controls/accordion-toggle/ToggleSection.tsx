import { Accordion, AccordionDetails, AccordionSummary, Checkbox, FormControlLabel, Radio } from '@mui/material';
import { FC, createContext, useCallback, useContext } from 'react';
import { useRecoilState } from 'recoil';

import { RecoilStateFamily } from '@/lib/recoil/types';

function useHandleCheckbox(onChecked: (checked: boolean) => void) {
  return useCallback((e, checked: boolean) => onChecked(checked), [onChecked]);
}

export const ToggleStateContext = createContext<RecoilStateFamily<boolean, string>>(null);

export const ToggleSectionGroup: FC<{ toggleState: RecoilStateFamily<boolean, string> }> = ({
  children,
  toggleState,
}) => {
  return <ToggleStateContext.Provider value={toggleState}>{children}</ToggleStateContext.Provider>;
};

interface ToggleSectionProps {
  id: string;
  label: string;
  forceSingle?: boolean;
  disabled?: boolean;
}

export const ToggleSection: FC<ToggleSectionProps> = ({
  id,
  label,
  forceSingle = false,
  disabled = false,
  children,
}) => {
  const toggleState = useContext(ToggleStateContext);
  const [show, setShow] = useRecoilState(toggleState(id));
  const handleShow = useHandleCheckbox(setShow);

  return (
    <Accordion disableGutters disabled={disabled} expanded={show} onChange={handleShow}>
      <AccordionSummary>
        <FormControlLabel
          control={
            forceSingle ? (
              <Radio checked={show} onChange={handleShow} />
            ) : (
              <Checkbox checked={show} onChange={handleShow} />
            )
          }
          label={label}
          onClick={
            // clicking on checkbox label shouldn't also trigger accordion change because then nothing happens
            (e) => e.preventDefault()
          }
        />
      </AccordionSummary>
      <AccordionDetails style={{ display: 'block' }}>{children}</AccordionDetails>
    </Accordion>
  );
};
