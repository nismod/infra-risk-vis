export const HAZARD_SOURCE = {
  getDataUrl({ hazardType, hazardParams: { returnPeriod, rcp, epoch, confidence } }, { scheme, range }) {
    const sanitisedRcp = rcp.replace('.', 'x');

    return `/raster/singleband/${hazardType}/${returnPeriod}/${sanitisedRcp}/${epoch}/${confidence}/{z}/{x}/{y}.png?colormap=${scheme}&stretch_range=[${range[0]},${range[1]}]`;
  },
};
