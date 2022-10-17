import { useCallback, useState } from 'react';

export function useTrigger() {
  const [trigger, setTrigger] = useState(0);

  const doTrigger = useCallback(() => {
    setTrigger((trigger) => trigger + 1);
  }, []);

  return [trigger, doTrigger] as const;
}
