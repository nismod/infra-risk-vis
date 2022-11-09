import { createContext, useContext } from 'react';

export const DataGroupContext = createContext<string>(null);

export function useDataGroup() {
  return useContext(DataGroupContext);
}

export function DataGroup({ group, children }) {
  return <DataGroupContext.Provider value={group}>{children}</DataGroupContext.Provider>;
}
