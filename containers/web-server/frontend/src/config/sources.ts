function paramString(params: string[]) {
  return params.length ? `?${params.join('&')}` : '';
}

export const SOURCES = {
  raster: {
    getUrl: ({ path, scheme, range }) => {
      const params = [];
      if (scheme) {
        params.push(`colormap=${scheme}`);
      }
      if (range) {
        params.push(`stretch_range=[${range[0]},${range[1]}]`);
      }
      return `/api/tiles/${path}/{z}/{x}/{y}.png${paramString(params)}`;
    },
  },
  vector: {
    getUrl: (datasetId) => `/vector/data/${datasetId}.json`,
  },
};

export type SourceName = keyof typeof SOURCES;
