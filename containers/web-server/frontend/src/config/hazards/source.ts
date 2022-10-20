import { SOURCES } from '../sources';

export const HAZARD_SOURCE = {
  getDataUrl({ hazardType, hazardParams }, { scheme, range }) {
    let path: string;
    // TODO: Gather required keys from API
    if (hazardType === 'earthquake') {
      const { returnPeriod, medium } = hazardParams;
      path = `earthquake/${returnPeriod}/${medium}`;
    } else {
      const { returnPeriod, rcp, epoch, gcm } = hazardParams;
      const sanitisedRcp = rcp.replace('.', 'x');

      if (hazardType === 'cyclone') {
        path = `${hazardType}/${returnPeriod}/${gcm}`;
      } else if (hazardType === 'extreme_heat_occurrence') {
        // TODO: Add support for exposure metric (as well as occurrence)
        path = `extreme_heat/occurrence/${sanitisedRcp}/${epoch}/${gcm}`;
      } else if (hazardType === 'extreme_heat_exposure') {
        // TODO: Add support for exposure metric (as well as occurrence)
        path = `extreme_heat/exposure/${sanitisedRcp}/${epoch}/${gcm}`;
      } else {
        path = `${hazardType}/${returnPeriod}/${sanitisedRcp}/${epoch}/${gcm}`;
      }
    }

    return SOURCES.raster.getUrl({
      path,
      scheme,
      range,
    });
  },
};
