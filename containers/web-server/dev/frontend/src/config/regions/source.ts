export const REGIONS_SOURCE = {
  getDataUrl({ regionLevel }) {
    return `/vector/data/regions_${regionLevel}.json`;
  },
};
