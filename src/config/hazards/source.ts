export const HAZARD_SOURCE = {
  getDataUrl({ hazardType, hazardParams: { returnPeriod, rcp, epoch, gcm } }, { scheme, range }) {
    const sanitisedRcp = rcp.replace('.', 'x');

    return `/raster/singleband/${hazardType}/${returnPeriod}/${sanitisedRcp}/${epoch}/${gcm}/{z}/{x}/{y}.png?colormap=${scheme}&stretch_range=[${range[0]},${range[1]}]`;
  },
};
