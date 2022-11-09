import { useRecoilValue } from 'recoil';

import { paramOptionsState, paramValueState, paramsState, useUpdateDataParam } from '@/state/data-params';

export const DataParam = ({ group, id, children }) => {
  const config = useRecoilValue(paramsState(group));

  const value = useRecoilValue(paramValueState({ group, param: id }));
  const updateValue = useUpdateDataParam(group, id);
  const options = useRecoilValue(paramOptionsState({ group, param: id }));

  if (config == null) return null;

  return typeof children === 'function' ? children({ value: value, onChange: updateValue, options }) : children;
};
