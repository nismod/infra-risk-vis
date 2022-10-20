export const SOURCES = {
  raster: {
    getUrl: ({ path, scheme, range }) =>
      `/api/tiles/${path}/{z}/{x}/{y}.png?colormap=${scheme}&stretch_range=[${range[0]},${range[1]}]`,
  },
};

export type SourceName = keyof typeof SOURCES;
