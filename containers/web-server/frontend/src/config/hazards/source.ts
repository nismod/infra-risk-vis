export const HAZARD_SOURCE = {
  getDataUrl({ hazardType, hazardParams: { returnPeriod, rcp, epoch, gcm } }, { scheme, range }) {
    const sanitisedRcp = rcp.replace('.', 'x');
    // TODO: Swap this to parse the meta from API about DB source instead of hardcoded
    return `/raster/singleband/aqueduct/${hazardType}/${returnPeriod}/${sanitisedRcp}/${epoch}/${gcm}/{z}/{x}/{y}.png?colormap=${scheme}&stretch_range=[${range[0]},${range[1]}]`;
  },
};
