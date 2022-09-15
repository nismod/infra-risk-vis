export type ConfigTree<T> = (ConfigTree<T> | false | null | undefined | T)[];
