import { useRecoilValue } from 'recoil';
import { dataParamOptionsState, dataParamState, useUpdateDataParam } from 'state/data-params';

export const DataParam = ({ group, id, children }) => {
  const value = useRecoilValue(dataParamState({ group, param: id }));
  const updateValue = useUpdateDataParam(group, id);
  const options = useRecoilValue(dataParamOptionsState({ group, param: id }));

  return typeof children === 'function' ? children({ value: value, onChange: updateValue, options }) : children;
};
