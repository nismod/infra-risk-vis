import { makeConfig } from '../helpers';

export const SOURCES = makeConfig([
  // {
  //   id: 'tileservergl',
  //   getUrl: ({ dataset }) => `http://localhost:8080/data/${dataset}.json`,
  // },
  // {
  //   id: 'terracotta',
  //   getUrl: ({ dataset }) => `http://localhost:8080/data/${dataset}.json`,
  // },
]);

export type SourceName = keyof typeof SOURCES;
