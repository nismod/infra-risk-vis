import { Paper, Table, TableBody, TableCell, TableContainer, TableHead, TableRow } from '@mui/material';
import React from 'react';

import ScrollToTop from '@/lib/react/ScrollToTop';

export const DataPage = () => (
  <article>
    <ScrollToTop />
    <h1>Data Sources and Access</h1>

    <h4>Disclaimer</h4>
    <p>
      This website was created for communication purposes. The information included here must not be used for the design
      of hazard-resistant structures or to support any important decision involving human life, capital or property. The
      values of risk or hazard in this map do not constitute an alternative nor do they replace building actions defined
      in national building codes. Readers seeking this information should consult national databases. The data presented
      here are the combination of results computed using multiple input models covering the majority of landmass. These
      models represent the best information currently publicly accessible, as far as the authors are aware. The
      presentation here results from an integration process that is solely the responsibility of the University of
      Oxford.
    </p>

    <h2>Hazard Data</h2>

    <TableContainer component={Paper}>
      <Table aria-label="simple table">
        <TableHead>
          <TableRow>
            <TableCell>Dataset</TableCell>
            <TableCell>Source</TableCell>
            <TableCell>Citation</TableCell>
            <TableCell>License</TableCell>
            <TableCell>Notes</TableCell>
          </TableRow>
        </TableHead>
        <TableBody>
          <TableRow>
            <TableCell>Coastal and River flooding</TableCell>
            <TableCell>
              <a href="https://www.wri.org/data/aqueduct-floods-hazard-maps">
                WRI Aqueduct Floods Hazard Maps
              </a>
            </TableCell>
            <TableCell>

              Ward, P.J., H.C. Winsemius, S. Kuzma, M.F.P. Bierkens, A. Bouwman,
              H. de Moel, A. Díaz Loaiza, et al. 2020. “Aqueduct Floods
              Methodology.” Technical Note. Washington, D.C.: World Resources
              Institute. Available online at: <a
              href="https://www.wri.org/publication/aqueduct-floods-methodology">
              www.wri.org/publication/aqueduct-floods-methodology</a>.

            </TableCell>
            <TableCell>

              All the products, methodologies, and datasets that make up
              Aqueduct are available for use under the{' '} <a
              href="https://creativecommons.org/licenses/by/4.0/"> Creative
              Commons Attribution International 4.0 License </a>.

            </TableCell>
            <TableCell>

              Inundation depth in meters for coastal and riverine floods over
              1km grid squares. 1 in 2 to 1 in 1000 year return periods.
              Baseline, RCP 4.5 &amp; 8.5 emission scenarios. Current and future
              maps in 2030, 2050 and 2080.

            </TableCell>
          </TableRow>
          <TableRow>
            <TableCell>Extreme Heat and Drought</TableCell>
            <TableCell>
              <a href="https://data.isimip.org/search/tree/ISIMIP2b/DerivedOutputData/lange2020/">
                Lange et al 2020, ISIMIP
              </a>
            </TableCell>
            <TableCell>

              Lange, S., Volkholz, J., Geiger, T., Zhao, F., Vega, I., Veldkamp,
              T., et al. (2020). Projecting exposure to extreme climate impact
              events across six event categories and three spatial scales.
              Earth's Future, 8, e2020EF001616. <a
              href="https://doi.org/10.1029/2020EF001616">DOI
              10.1029/2020EF001616</a>

            </TableCell>
            <TableCell>CC0 1.0</TableCell>
            <TableCell>

              Annual probability of drought (soil moisture below a baseline
              threshold) or extreme heat (temperature and humidity-based
              indicators over a threshold) events on a 0.5° grid. 8 hydrological
              models forced by 4 GCMs under baseline, RCP 2.6 &amp; 6.0 emission
              scenarios. Current and future maps in 2030, 2050 and 2080. <br/>
              <br/>

              The ISIMIP2b climate input data and impact model output data
              analyzed in this study are available in the ISIMIP data repository
              at ESGF, see
              https://esg.pik-potsdam.de/search/isimip/?project=ISIMIP2b&product=input
              and
              https://esg.pik-potsdam.de/search/isimip/?project=ISIMIP2b&product=output,
              respectively. More information about the GHM, GGCM, and GVM output
              data is provided by Gosling et al. (2020), Arneth et al. (2020),
              and Reyer et al. (2019), respectively. <br/> <br/>

              Event definitions are given in Lange et al, table 1. Land area is
              exposed to drought if monthly soil moisture falls below the 2.5th
              percentile of the preindustrial baseline distribution for at least
              seven consecutive months. Land area is exposed to extreme heat if
              both a relative indicator based on temperature (Russo et al 2015,
              2017) and an absolute indicator based on temperature and relative
              humidity (Masterton &amp; Richardson, 1979) exceed their
              respective threshold value. <br/> <br/>

              Note that the time series of extreme events given by Lange et al
              has been processed into an annual probability of occurrence by
              the GRI team for visualisation purposes.

            </TableCell>
          </TableRow>
          <TableRow>
            <TableCell>Tropical Cyclones</TableCell>
            <TableCell>

                STORM Tropical Cyclone Maximum Windspeeds, <a
                href="https://data.4tu.nl/articles/dataset/STORM_tropical_cyclone_wind_speed_return_periods/12705164/3">
                Present </a>{" and "} <a
                href="https://data.4tu.nl/articles/dataset/STORM_climate_change_tropical_cyclone_wind_speed_return_periods/14510817/3">
                Future climate </a>.

            </TableCell>
            <TableCell>

              Bloemendaal, Nadia; de Moel, H. (Hans); Muis, S; Haigh, I.D.
              (Ivan); Aerts, J.C.J.H. (Jeroen) (2020): STORM tropical cyclone
              wind speed return periods. 4TU.ResearchData. Dataset. <a
              href="https://doi.org/10.4121/12705164.v3">DOI
              10.4121/12705164.v3</a> and Bloemendaal, Nadia; de Moel, Hans;
              Dullaart, Job; Haarsma, R.J. (Reindert); Haigh, I.D. (Ivan);
              Martinez, Andrew B.; et al. (2022): STORM climate change tropical
              cyclone wind speed return periods. 4TU.ResearchData. Dataset. <a
              href="https://doi.org/10.4121/14510817.v3">DOI
              10.4121/14510817.v3</a>

            </TableCell>
            <TableCell>CC0 1.0</TableCell>
            <TableCell>

              Tropical cyclone maximum wind speed (in m/s) return periods,
              generated using the STORM climate change datasets. 1 in 10 to 1 in
              10,000 year return periods at 10 km resolution. Baseline and RCP
              8.5 climate scenarios. Current and future maps in 2030, 2050 and
              2080.

            </TableCell>
          </TableRow>
          <TableRow>
            <TableCell>Seismic Risk</TableCell>
            <TableCell>
              <a href="https://www.globalquakemodel.org/gem-maps/global-earthquake-hazard-map">
                GEM Global Earthquake Hazard Map
              </a>
            </TableCell>
            <TableCell>

              Pagani M, Garcia-Pelaez J, Gee R, Johnson K, Silva V, Simionato M,
              Styron R, Vigano D, Danciu L, Monelli D, Poggi V, Weatherill G.
              (2019). The 2018 version of the Global Earthquake Model: Hazard
              component. Earthquake Spectra, 36(1), <a
              href="https://doi.org/10.1177/8755293020931866">DOI:
              10.1177/8755293020931866</a>.

            </TableCell>
            <TableCell>

              The visualisation of seismic hazard data is licensed under the
              terms of the{' '} <a
              href="https://creativecommons.org/licenses/by-nc-sa/4.0/">
              Creative Commons Attribution-NonCommercial-ShareAlike 4.0
              International License (CC BY-NC-SA) </a> . The underlying data are
              confidential,{' '} <a
              href="mailto:info@globalquakemodel.org">contact GEM</a>{' '} for
              more information.

            </TableCell>
            <TableCell>

              The Global Earthquake Model (GEM) Global Seismic Hazard Map
              (version 2018.1) depicts the geographic distribution of the Peak
              Ground Acceleration (PGA) with a 10% probability of being exceeded
              in 50 years, computed for reference rock conditions (shear wave
              velocity, VS30, of 760-800 m/s).

            </TableCell>
          </TableRow>

        </TableBody>
      </Table>
    </TableContainer>

    <h2>Exposure Data</h2>

    <TableContainer component={Paper}>
      <Table aria-label="simple table">
        <TableHead>
          <TableRow>
            <TableCell>Dataset</TableCell>
            <TableCell>Source</TableCell>
            <TableCell>Citation</TableCell>
            <TableCell>License</TableCell>
            <TableCell>Notes</TableCell>
          </TableRow>
        </TableHead>
        <TableBody>

          <TableRow>
            <TableCell>

               Roads and Rail

            </TableCell>
            <TableCell>

              <a href="https://planet.openstreetmap.org/">OpenStreetMap</a>

            </TableCell>
            <TableCell>

              <a href="https://www.openstreetmap.org/copyright">© OpenStreetMap contributors https://www.openstreetmap.org/copyright</a>

            </TableCell>
            <TableCell>

              ODbL

            </TableCell>
            <TableCell>

              Extract from OpenStreetMap October 2021. All roads tagged as
              trunk, motorway, primary, secondary or tertiary, all rail lines
              tagged as rail and railway stations.

            </TableCell>
          </TableRow>

          <TableRow>
            <TableCell>

              <a href="https://doi.org/10.5281/zenodo.3628142">Gridfinder Power
              Transmission lines</a>

            </TableCell>
            <TableCell>

              Source

            </TableCell>
            <TableCell>

              Arderne, C., Zorn, C., Nicolas, C. et al. Predictive mapping of
              the global power system using open data. Sci Data 7, 19 (2020).
              https://doi.org/10.1038/s41597-019-0347-4


            </TableCell>
            <TableCell>

              CC BY 4.0

            </TableCell>
            <TableCell>

              Predicted distribution and transmission line network, with
              existing OpenStreetMap lines tagged in the 'source' column and ©
              OpenStreetMap contributors


            </TableCell>
          </TableRow>


          <TableRow>
            <TableCell>

              Power plants

            </TableCell>
            <TableCell>

              <a
              href="https://datasets.wri.org/dataset/globalpowerplantdatabase">WRI
              Global Powerplants Database</a>

            </TableCell>
            <TableCell>

              Global Energy Observatory, Google, KTH Royal Institute of
              Technology in Stockholm, Enipedia, World Resources Institute.
              2018. Global Power Plant Database. Published on Resource Watch and
              Google Earth Engine; http://resourcewatch.org/
              https://earthengine.google.com/

            </TableCell>
            <TableCell>

              CC BY 4.0

            </TableCell>
            <TableCell>

              The Global Power Plant Database is a comprehensive, open source
              database of power plants around the world. It centralizes power
              plant data to make it easier to navigate, compare and draw
              insights for one's own analysis. The database covers approximately
              35,000 power plants from 167 countries and includes thermal plants
              (e.g. coal, gas, oil, nuclear, biomass, waste, geothermal) and
              renewables (e.g. hydro, wind, solar). Each power plant is
              geolocated and entries contain information on plant capacity,
              generation, ownership, and fuel type.

            </TableCell>
          </TableRow>


          <TableRow>
            <TableCell>

              Population and built-up area

            </TableCell>
            <TableCell>

              JRC Global Human Settlement Layer

            </TableCell>
            <TableCell>

              Dataset: Pesaresi M., Politis P. (2022): GHS built-up surface
              grid, derived from Sentinel2 composite and Landsat, multitemporal
              (1975-2030)European Commission, Joint Research Centre (JRC) PID:
              http://data.europa.eu/89h/d07d81b4-7680-4d28-b896-583745c27085,
              doi:10.2905/D07D81B4-7680-4D28-B896-583745C27085<br /><br/>

              Concept & Methodology: Schiavina M., Melchiorri M., Pesaresi M.,
              Politis P., Freire S., Maffenini L., Florio P., Ehrlich D., Goch
              K., Tommasi P., Kemper T. GHSL Data Package 2022, JRC 129516, ISBN
              978-92-76-53071-8 doi:10.2760/19817<br /><br/>

              Dataset: Schiavina M., Freire S., MacManus K. (2022): GHS-POP
              R2022A - GHS population grid multitemporal (1975-2030).European
              Commission, Joint Research Centre (JRC) PID:
              http://data.europa.eu/89h/d6d86a90-4351-4508-99c1-cb074b022c4a,
              doi:10.2905/D6D86A90-4351-4508-99C1-CB074B022C4A<br /><br/>

              Concept & Methodology: Freire S., MacManus K., Pesaresi M.,
              Doxsey-Whitfield E., Mills J. (2016) Development of new open and
              free multi-temporal global population grids at 250 m resolution.
              Geospatial Data in a Changing World; Association of Geographic
              Information Laboratories in Europe (AGILE), AGILE 2016


            </TableCell>
            <TableCell>

              <a href="https://ec.europa.eu/info/legal-notice_en#copyright-notice">CC-BY 4.0</a>

            </TableCell>
            <TableCell>

              GHS-POP R2022A - GHS population grid multitemporal (1975-2030),
              epoch: 2020, resolution: 1 km, coordinate system: Mollweide,
              reprojected for visualisation. For GHS-POP (GHS-POP_GLOBE_R2022A),
              the Sentinel/Landsat based GHS-BUILT-S (GHS-BUILT-S_GLOBE_R2022A,
              version 1.0) was used as target for disaggregation of population
              estimates. The base source for population estimates (both census
              unit counts and geometries) was the raw dataset (census population
              at the census year and growth rates) of the Gridded Population of
              the World, version 4.11 (GPWv4.11), from CIESIN/SEDAC
              (https://sedac.ciesin.columbia.edu/data/collection/gpw-v4/whatsnewrev11).<br /><br/>

              GHS-BUILT-S R2022A - GHS built-up surface grid, derived from
              Sentinel-2 composite and Landsat, multitemporal (1975-2030),
              epoch: 2020, resolution: 1 km, coordinate system: Mollweide,
              classification: Total RES+NRES or Non Residential, reprojected for
              visualisation. The GHS-BUILT-S spatial raster dataset depicts the
              distribution of the built-up (BU) surfaces estimates between 1975
              and 2030 in 5 years intervals and two functional use components a)
              the total BU surface and b) the non-residential (NRES) BU surface.
              The data is made by spatial-temporal interpolation of five
              observed collections of multiple-sensor, multiple-platform
              satellite imageries: Landsat (MSS, TM, ETM sensor) supports the
              1975, 1990, 2000, and 2014 epochs, while Sentinel-2 (S2) composite
              (GHS-composite-S2 R2020A) supports the 2018 epoch.

          </TableCell>
          </TableRow>

          <TableRow>
            <TableCell>

              Health site locations


            </TableCell>
            <TableCell>

              <a href="http://healthsites.io">healthsites.io</a>

            </TableCell>
            <TableCell>

              This data was generated as an extract from the OpenStreetMap global
              open database (http://openstreetmap.org) by the Healthsites.io
              (http://healthsites.io) project. Data: Open Database License
              http://opendatacommons.org/licenses/odbl/ Data credits: ©
              OpenStreetMap contributors http://www.openstreetmap.org/copyright

            </TableCell>
            <TableCell>

              ODbL

            </TableCell>
            <TableCell>

              Health site locations as points/polygons

            </TableCell>
          </TableRow>


          <TableRow>
            <TableCell>

               Cement and Steel Production Assets

            </TableCell>
            <TableCell>

              <a
              href="https://www.cgfi.ac.uk/spatial-finance-initiative/database-downloads/">Global
              Databases of Cement and Iron and Steel Production Assets, Spatial
              Finance Initiative</a>

            </TableCell>
            <TableCell>

              McCarten, M., Bayaraa, M., Caldecott, B., Christiaen, C., Foster,
              P., Hickey, C., Kampmann, D., Layman, C., Rossi, C., Scott, K.,
              Tang, K., Tkachenko, N., and Yoken, D. 2021. Global Database of
              Cement Production Assets. Spatial Finance Initiative.<br/><br/>

              McCarten, M., Bayaraa, M., Caldecott, B., Christiaen, C., Foster,
              P., Hickey, C., Kampmann, D., Layman, C., Rossi, C., Scott, K.,
              Tang, K., Tkachenko, N., and Yoken, D., 2021. Global Database of
              Iron and Steel Production Assets. Spatial Finance Initiative


            </TableCell>
            <TableCell>

              CC BY 4.0

            </TableCell>
            <TableCell>

              Cement and Steel Asset site locations as points

            </TableCell>
          </TableRow>

          <TableRow>
            <TableCell>

              Soil Organic Carbon stock

            </TableCell>
            <TableCell>

              <a href="https://soilgrids.org/">SoilGrids 2.0</a>

            </TableCell>
            <TableCell>

              Poggio, L., de Sousa, L.M., Batjes, N.H., Heuvelink, G.B.M.,
              Kempen, B., Ribeiro, E., Rossiter, D., 2021. SoilGrids 2.0:
              producing soil information for the globe with quantified spatial
              uncertainty. SOIL 7, 217–240.
              https://doi.org/10.5194/soil-7-217-2021


            </TableCell>
            <TableCell>

              CC-BY 4.0

            </TableCell>
            <TableCell>


              Soil organic carbon content at 0-30cm, in tonnes/hectare,
              aggregated to 1000m grid.

              Soil organic carbon content (fine earth fraction)
              in dg/kg at 6 standard depths. Predictions were derived using a
              digital soil mapping approach based on Quantile Random Forest,
              drawing on a global compilation of soil profile data and
              environmental layers. This map is the result of resampling the
              mean SoilGrids 250 m predictions (Poggio et al. 2021) for each
              1000 m cell.

            </TableCell>
          </TableRow>

        </TableBody>
      </Table>
    </TableContainer>

    <h2>Vulnerability Data</h2>

    <TableContainer component={Paper}>
      <Table aria-label="simple table">
        <TableHead>
          <TableRow>
            <TableCell>Dataset</TableCell>
            <TableCell>Source</TableCell>
            <TableCell>Citation</TableCell>
            <TableCell>License</TableCell>
            <TableCell>Notes</TableCell>
          </TableRow>
        </TableHead>
        <TableBody>

          <TableRow>
            <TableCell>Access to Healthcare</TableCell>
            <TableCell>

              <a href="https://www.nature.com/articles/s41591-020-1059-1">Global
              maps of travel time to healthcare facilities</a>

            </TableCell>
            <TableCell>

              Weiss, D.J., Nelson, A., Vargas-Ruiz, C.A. et al. Global maps of
              travel time to healthcare facilities. Nat Med 26, 1835–1838
              (2020). https://doi.org/10.1038/s41591-020-1059-1

            </TableCell>
            <TableCell>CC BY 4.0</TableCell>
            <TableCell>

              Motorised/non-motorised travel time on 30arcsec grid.

            </TableCell>
          </TableRow>

          <TableRow>
            <TableCell>Human Development</TableCell>
            <TableCell>

              <a href="https://globaldatalab.org/shdi/">Global Data Lab
              Sub-national human development indices</a>

            </TableCell>
            <TableCell>

              Global Data Lab (2019) Subnational Human Development Index (SHDI)
              Available at https://globaldatalab.org/shdi/. The SHDI is an
              average of the subnational values of three dimensions: education,
              health and standard of living. To compute the SHDI on the basis of
              the three dimension indices, the geometric mean of the three
              indices is taken. Three major data sources were used to create the
              SHDI database: statistical offices (including Eurostat, the
              statistical office of the European Union), the Area Database of
              the Global Data Lab, and data from the HDI website of the Human
              Development Report Office of the United Nations Development
              Program. Given that household surveys and censuses are not held
              every year, for many countries the indicators are only available
              for a restricted number of years. To obtain their values for the
              whole period 1990–2017, the missing information was estimated by
              interpolation or extrapolation techniques. This estimation process
              was facilitated by the fact that the UNDP Database contains the
              national values for all four indicators for each year in this
              period, which means that only the subnational variation had to be
              interpolated or extrapolated. For a complete list of sources and
              surveys used, please refer to the Area Database's Data Sources
              page.

            </TableCell>
            <TableCell>Free for use with acknowledgement of data source https://globaldatalab.org/termsofuse/</TableCell>
            <TableCell>

              Development, Health, Education and Income indices for 186
              countries, 1783 sub-national regions

            </TableCell>
          </TableRow>
          <TableRow>
            <TableCell>Biodiversity</TableCell>
            <TableCell>

              <a
              href="https://data.nhm.ac.uk/dataset/global-map-of-the-biodiversity-intactness-index-from-newbold-et-al-2016-science/resource/8531b4dc-bd44-4586-8216-47b3b8d60e85">Biodiversity
              Intactness Index</a>

            </TableCell>
            <TableCell>

              Tim Newbold; Lawrence Hudson; Andy Arnell; Sara Contu et al.
              (2016). Map of Biodiversity Intactness Index (from Global map of
              the Biodiversity Intactness Index, from Newbold et al. (2016)
              Science) [Data set resource]. Natural History Museum.
              https://data.nhm.ac.uk/dataset/global-map-of-the-biodiversity-intactness-index-from-newbold-et-al-2016-science/resource/8531b4dc-bd44-4586-8216-47b3b8d60e85

            </TableCell>
            <TableCell>CC BY 4.0</TableCell>
            <TableCell>

              3 arcsec grid

            </TableCell>
          </TableRow>
          <TableRow>
            <TableCell>Protected Areas</TableCell>
            <TableCell>

              <a
              href="https://www.protectedplanet.net/en/thematic-areas/wdpa?tab=WDPA">World
              Database of Protected Areas</a>

            </TableCell>
            <TableCell>

              UNEP-WCMC and IUCN (2022), Protected Planet: The World Database on
              Protected Areas (WDPA) [On-line], [October 2022], Cambridge, UK:
              UNEP-WCMC and IUCN. Available at: www.protectedplanet.net.

            </TableCell>
            <TableCell>No Commercial Use, No Reposting and/or Redistribution without written consent</TableCell>
            <TableCell>

              Protected area locations as points/polygons

            </TableCell>
          </TableRow>
          <TableRow>
            <TableCell>Forest Integrity</TableCell>
            <TableCell>

              <a href="https://www.nature.com/articles/s41467-020-19493-3">Forest Landscape Integrity Index</a>

            </TableCell>
            <TableCell>

              Grantham, H.S., Duncan, A., Evans, T.D. et al. Anthropogenic modification of forests means only 40% of remaining forests have high ecosystem integrity. Nat Commun 11, 5978 (2020). https://doi.org/10.1038/s41467-020-19493-3 Data are available at www.forestlandscapeintegrity.com. The datasets used to develop the Forest Landscape Integrity Index can be found at the following websites: tree cover and loss http://earthenginepartners.appspot.com/science-2013-global-forest, tree cover loss driver https://data.globalforestwatch.org/datasets/f2b7de1bdde04f7a9034ecb363d71f0e, potential forest cover https://data.globalforestwatch.org/datasets/potential-forest-coverage ESA-CCI Land Cover https://maps.elie.ucl.ac.be/CCI/viewer/index.php Open Street Maps https://www.openstreetmap.org, croplands https://lpdaac.usgs.gov/news/release-of-gfsad-30-meter-cropland-extent-products/, surface water https://global-surface-water.appspot.com/, protected areas https://www.protectedplanet.net/en.

            </TableCell>
            <TableCell>Published as available with article (license not specified)</TableCell>
            <TableCell>

              10 arcsec grid

            </TableCell>
          </TableRow>
        </TableBody>
      </Table>
    </TableContainer>


    <h2>Contextual Map Data</h2>

    <p>
      Background map data is &copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors,
      style &copy; <a href="https://carto.com/attributions">CARTO</a>.
    </p>

    <p>
      Satellite imagery background is derived from{' '}
      <a href="https://s2maps.eu">Sentinel-2 cloudless - https://s2maps.eu</a> by{' '}
      <a href="https://eox.at">EOX IT Services GmbH</a> (Contains modified Copernicus Sentinel data 2020).
    </p>

    <p>
      Photo credit: Hurricane Irma, 7 September 2017. Data: MODIS/Terra (NASA
      WorldView). Processed by Antti Lipponen (<a
      href="https://twitter.com/anttilip">@anttilip</a>) <a
      href="https://creativecommons.org/licenses/by/2.0/">CC-BY</a>
    </p>


  </article>
);
