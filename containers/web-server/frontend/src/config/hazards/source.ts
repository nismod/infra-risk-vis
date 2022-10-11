export const HAZARD_SOURCE = {
  getDataUrl({ hazardType, hazardParams: { returnPeriod, rcp, epoch, gcm } }, { scheme, range }) {
    const sanitisedRcp = rcp.replace('.', 'x');
    // TODO: Gather required keys from API
    if (hazardType === 'cyclone') {
      return `/api/tiles/${hazardType}/${returnPeriod}/${gcm}//{z}/{x}/{y}.png?colormap=${scheme}&stretch_range=[${range[0]},${range[1]}]`;
    }
    else if (hazardType === 'extreme_heat') {
      // TODO: Add support for exposure metric (as well as occurrence)
      return `/api/tiles/${hazardType}/occurrence/${sanitisedRcp}/${epoch}/${gcm}//{z}/{x}/{y}.png?colormap=${scheme}&stretch_range=[${range[0]},${range[1]}]`;
    } else {
      return `/api/tiles/${hazardType}/${returnPeriod}/${sanitisedRcp}/${epoch}/${gcm}/{z}/{x}/{y}.png?colormap=${scheme}&stretch_range=[${range[0]},${range[1]}]`;
    }
  },
};
