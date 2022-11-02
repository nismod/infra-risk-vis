import { SOURCES } from '../sources';

export const HAZARD_SOURCE = {
  getDataUrl({ hazardType, hazardParams }, { scheme, range }) {
    let path: string;
    // TODO: Gather required keys from API
    if (hazardType === 'earthquake') {
      const { rp, medium } = hazardParams;
      path = `earthquake/${rp}/${medium}`;
    } else {
      const { rp, rcp, epoch, gcm } = hazardParams;
      const sanitisedRcp = rcp?.replace('.', 'x');

      if (hazardType === 'cyclone') {
        path = `${hazardType}/${rp}/${gcm}`;
      } else if (hazardType === 'extreme_heat') {
        // TODO: Add support for exposure metric (as well as occurrence)
        path = `extreme_heat/occurrence/${sanitisedRcp}/${epoch}/${gcm}`;
        /*
        else if (hazardType === 'extreme_heat_exposure') {
          path = `extreme_heat/exposure/${sanitisedRcp}/${epoch}/${gcm}`;
        }
        */
      } else if (hazardType === 'drought') {
        path = `drought/occurrence/${sanitisedRcp}/${epoch}/${gcm}`;
      } else {
        path = `${hazardType}/${rp}/${sanitisedRcp}/${epoch}/${gcm}`;
      }
    }

    return SOURCES.raster.getUrl({
      path,
      scheme,
      range,
    });
  },
};
