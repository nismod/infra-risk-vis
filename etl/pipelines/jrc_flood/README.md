# JRC Global River Flood

Global river flood hazard maps is a gridded data set representing inundation
along the river network, for seven different flood return periods. The input
river flow data for the new maps are produced by means of the open-source
hydrological model LISFLOOD, while inundation simulations are performed with the
hydrodynamic model LISFLOOD-FP.

## Citation

Baugh, Calum; Colonese, Juan; D'Angelo, Claudia; Dottori, Francesco; Neal,
Jeffrey; Prudhomme, Christel; Salamon, Peter (2024): Global river flood hazard
maps. European Commission, Joint Research Centre (JRC) [Dataset] PID:
http://data.europa.eu/89h/jrc-floods-floodmapgl_rp50y-tif

## Metadata

GENERAL INFORMATION

G01. Dataset version 2.0.0
G02. Name of the dataset: Global river flood hazard maps
G03. Description of the dataset: Global river flood hazard maps is a gridded data set representing inundation along the river network, for seven different flood return periods. The input river flow data for the new maps are produced by means of the open-source hydrological model LISFLOOD, while inundation simulations are performed with the hydrodynamic model LISFLOOD-FP.
G04. Creator of data set: Copernicus Emergency Management Service
G05. DOI of dataset: NA
G06. PID of dataset: NA
G07. Keywords (author defined): riverine flood, hazard, global
G08. Language(s) used in the dataset: English
G09. Last update of the README file: 14.03.2024

DATA DESCRIPTION

D01. Horizontal coverage: Global
D02. Horizontal resolution: 3 arc seconds (~90 m)
D03. Spatial gaps: Only land areas are covered by this data set
D04. Temporal coverage: NA
D05. Temporal resolution: NA
D06. Temporal gaps: NA
D07. Number of available variables: 1 (flood water depth)
D08. Variables available at daily resolution: NA
D09. Variables available at 6-hourly resolution: NA
D10. Units: meters
D11. Update frequency: irregular
D12. Projection: wgs_1984
D13. Data type: 3 arc seconds (~90m) grid
D14. Available version(s): 2

PROJECT INFORMATION

P01. Project information: Global river flood hazard maps have been produced as part of the Global Flood Awareness System (GloFAS) of the Copernicus Emergency Management Service.
P02. Project website: https://emergency.copernicus.eu/
P03. Project funder: European Commission

FILE OVERVIEW

F01. Number of files described by the README-file: 266 files per return period. Files are stored in folders for each return period.
F02. The different files represent the tiles for the flood hazard maps for different return periods and the permanent water bodies used to patch the flood hazard maps for permanent water bodies
F03. Naming conventions for file names: <IDXXX_N/SXX_E/WXX_RPXXX>\_depth.tif
F04. Explanation of abbreviations: RP = Return period, XXX denotes the return period year or the tile identifier, N/S/W/E - North/South/East/West;
F05. File formats: Tagged Image File Format tif
F06. Software necessary to open the file: any software that opens tif, e.g. ArcGIS, QGIS etc.

STORAGE INFORMATION

S01. Where are the data stored?: JRC Data Store, https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/CEMS-GLOFAS/flood_hazard/

DATA ACCESS AND SHARING

A01. Recommended citation for the dataset: Baugh, Calum; Colonese, Juan; D'Angelo, Claudia; Francesco, Dottori; Neal, Jeffrey; Prudhomme; Christe; Salamon, Peter (2024): Modelled flood inundation for different return period scenarios at the global scale. European Commission, Joint Research Centre (JRC)
A02. License information, restrictions on use: no restrictions, free and open Copernicus product
A03. Copyright statement: https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/CEMS-GLOFAS/copyright.txt
A04. User support service: contact person outlined in the JRC Data Catalogue

RELATIONSHIPS

R01. Reference publication of this dataset: Baugh, Calum; Colonese, Juan; D'Angelo, Claudia; Francesco, Dottori; Neal, Jeffrey; Prudhomme, Christel; Salamon, Peter. (2024): An updated dataset of global and European flood hazard maps. Manuscript under preparation.
R02. Grimaldi, Stefania; Salamon, Peter; Disperati, Juliana; Zsoter, Ervin; Russo, Carlo; Ramos, Arthur; Carton, Corentin; Barnard, Chris; Hansford, Eleanor; Gomes, Goncalo; Prudhomme, Christel (2022): GloFAS v4.0 hydrological reanalysis. European Commission, Joint Research Centre (JRC) [Dataset] PID: http://data.europa.eu/89h/f96b7a19-0133-4105-a879-0536991ca9c5
R03. This dataset incorporates the following other datasets: MERIT-DEM, MERIT-HYDRO, GloFAS v4.0 hydrological reanalysis

## Copyright notice

(c) European Union, 1995-2024

The Commission's reuse policy is implemented by the Commission Decision of 12 December 2011
on the reuse of Commission documents [1]. Any copyright and/or sui generis right on the dataset
is licensed under the Creative Commons Attribution 4.0 International (CC BY 4.0) licence [2].
Reuse is allowed provided appropriate credit is given and any changes are indicated.

[1] https://eur-lex.europa.eu/eli/dec/2011/833/oj
[2] https://creativecommons.org/licenses/by/4.0

## Changelog

### V2.0.0 (2024-03-14)

Release of the updated global flood hazard maps. A description of the dataset
is currently being prepared as manuscript: Baugh, Calum; Colonese, Juan;
D'Angelo, Claudia; Francesco, Dottori; Neal, Jeffrey; Prudhomme, Christel;
Salamon, Peter. (2024): An updated dataset of global and European flood hazard
maps. Manuscript under preparation

### v1.0.0 (2016-05-15)

Initial release of the global flood hazard maps. A detailed description of the
maps and the methodology can be found in: Dottori F, Salamon P, Bianchi A,
Alfieri L, Hirpa F, Feyen L. Development and evaluation of a framework for
global flood hazard mapping. ADVANCES IN WATER RESOURCES 94; 2016. p. 87-102.
ISSN 0309-1708, https://doi.org/10.1016/j.advwatres.2016.05.002.
