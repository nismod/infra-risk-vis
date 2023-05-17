import { useEffect } from 'react';
import { atomFamily, useSetRecoilState } from 'recoil';

export const mobileTabHasContentState = atomFamily({
  key: 'mobileTabHasContentState',
  default: false,
});

/**
 * Use this component to indicate that a tab in the mobile UI version has content.
 * The `tabId` should match one of the `id` fields in `mobileTabsConfig` in this file.
 */
export const MobileTabContentWatcher = ({ tabId }) => {
  const setTabHasContent = useSetRecoilState(mobileTabHasContentState(tabId));

  useEffect(() => {
    setTabHasContent(true);

    return () => {
      setTabHasContent(false);
    };
  }, [setTabHasContent]);

  return null;
};
