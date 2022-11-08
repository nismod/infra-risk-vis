import { useCallback, useContext, useEffect } from 'react';
import { useRecoilCallback } from 'recoil';

import {
  VisibilityStateContext,
  usePathChildrenState,
  useVisibilityState,
} from '@/lib/data-selection/sidebar/components';
import { getSubPath, usePath } from '@/lib/data-selection/sidebar/paths';

const ChildVisibilityWatcher = ({ childPath, onVisibility }) => {
  const [visible] = useVisibilityState(childPath);

  useEffect(() => {
    onVisibility(visible);
  }, [onVisibility, visible]);

  return null;
};

export const EnforceSingleChild = () => {
  const path = usePath();
  const [subPaths] = usePathChildrenState(path);
  const visibilityState = useContext(VisibilityStateContext);

  const enforceLimit = useRecoilCallback(
    ({ set }) =>
      (newShown: string) => {
        for (const sp of subPaths) {
          if (sp !== newShown) {
            set(visibilityState(getSubPath(path, sp)), false);
          }
        }
      },
    [path, subPaths, visibilityState],
  );

  const handleChange = useCallback(
    (subPath, visibility) => {
      if (visibility) {
        enforceLimit(subPath);
      }
    },
    [enforceLimit],
  );

  return (
    <>
      {subPaths.map((sp) => (
        <ChildVisibilityWatcher key={sp} childPath={getSubPath(path, sp)} onVisibility={(v) => handleChange(sp, v)} />
      ))}
    </>
  );
};
