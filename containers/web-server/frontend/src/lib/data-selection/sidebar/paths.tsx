import { FC, createContext, useContext } from 'react';

export function getParentPath(path: string) {
  if (path === '') {
    throw new Error("Empty path doesn't have a parent");
  }
  const idx = path.lastIndexOf('/');
  if (idx === -1) {
    return '';
  }
  return path.slice(0, idx);
}

export function getSubPath(parentPath: string, path: string) {
  return parentPath === '' ? path : `${parentPath}/${path}`;
}

export const PathContext = createContext<string>('');

export function usePath() {
  return useContext(PathContext);
}

export const SubPath: FC<{ path: string }> = ({ path, children }) => {
  const parentPath = usePath();
  const subPath = getSubPath(parentPath, path);

  return <PathContext.Provider value={subPath}>{children}</PathContext.Provider>;
};
