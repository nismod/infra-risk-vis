import React from 'react';
import { Alert, Button, Paper, Table, TableBody, TableCell, TableContainer, TableHead, TableRow } from '@mui/material';
import ScrollToTop from 'lib/hooks/scroll-to-top';
import { ExtLink } from 'lib/nav';

export const GuidePage = () => (
  <article>
    <h1>User Guide for the Jamaica Systemic Risk Assessment Tool (J-SRAT)</h1>
    <p>
      <ExtLink href="/guide_media/J-SRAT User Guide_2022-05.pdf">Download PDF</ExtLink>
    </p>
    <h2 id="table-of-contents">Table of Contents</h2>
    <p>
      <a href="#how-to-use-this-guide">How to use this guide</a>
    </p>
    <p>
      <a href="#navigate-around-the-j-srat">Navigate around the J-SRAT</a>
    </p>
    <p>
      <a href="#explore-the-j-srat-data-and-analysis">Explore the J-SRAT data and analysis</a>
    </p>
    <ul>
      <li>
        <a href="#how-are-the-transport-energy-and-water-systems-represented">
          How are the transport, energy and water systems represented?
        </a>
      </li>
      <li>
        <a href="#how-are-flooding-and-hurricanes-represented">How are flooding and hurricanes represented?</a>
      </li>
      <li>
        <a href="#how-are-people-and-economic-activities-represented">
          How are people and economic activities represented?
        </a>
      </li>
      <li>
        <a href="#how-is-climate-risk-from-flooding-or-cyclones-represented">
          How is climate risk from flooding or cyclones represented?
        </a>
      </li>
      <li>
        <a href="#how-is-drought-risk-represented">How is drought risk represented?</a>
      </li>
      <li>
        <a href="#how-are-adaptation-options-presented">How are adaptation options presented?</a>
      </li>
      <li>
        <a href="#how-are-nature-based-solutions-represented">How are nature-based solutions represented?</a>
      </li>
    </ul>
    <p>
      <a href="#trace-the-climate-risk-and-adaptation-analysis-for-a-single-asset">
        Trace the climate risk and adaptation analysis for a single asset
      </a>
    </p>
    <p>
      <a href="#prioritise-and-evaluate-adaptation-interventions">Prioritise and evaluate adaptation interventions</a>
    </p>

    <p>
      <a href="#about-this-guide">About this guide</a>
    </p>

    <h2 id="how-to-use-this-guide">How to use this guide</h2>
    <a href="#table-of-contents">
      <small>Back to top</small>
    </a>
    <p>
      This guide is intended to accompany the interactive web-based platform, which is part of the Jamaica Systemic Risk
      Assessment Tool (J-SRAT), and visualises the results of a climate risk and adaptation analysis of Jamaica’s
      infrastructure.
    </p>
    <p>The overall objectives of the J-SRAT are to:</p>
    <ol type="1">
      <li>
        <p>
          Present the results of a climate risk analysis for Jamaica’s infrastructure networks (transport, energy,
          water) to estimate the economic impacts of physical climate risks and identify locations of vulnerability
          within infrastructure networks.
        </p>
      </li>
      <li>
        <p>
          Enable evaluation and prioritization of policies and investment options to reduce losses and enhance
          infrastructure resilience.
        </p>
      </li>
    </ol>
    <p>
      This guide aims to convey an initial understanding of the tool’s content and capabilities. First, we introduce how
      to navigate around the tool. Then we introduce each of the major data and results layers. Finally, we work through
      two more analytical use cases: stepping through the climate risk and adaptation analysis for a single asset, then
      prioritising adaptation interventions based on cost-benefit analysis and digging into the details to evaluate a
      particular intervention.
    </p>
    <p>
      For analysts, a detailed technical methodology report on the analysis is also available: Pant, R., Becher, O.,
      Haggis R., Majid, A., Russell, T., Verschuur, J., and Hall, J.W. (2022). Final technical report on methodology and
      implementation of the Jamaica Systemic Risk Assessment Tool (J-SRAT). Environmental Change Institute, University
      of Oxford, UK.
    </p>
    <p>
      For developers, the source code for the tool is developed and documented at{' '}
      <ExtLink href="https://github.com/nismod/infra-risk-vis/tree/release/jamaica">
        github.com/nismod/infra-risk-vis/tree/release/jamaica
      </ExtLink>
      . The analysis for Jamaica is produced using the code and models at{' '}
      <ExtLink href="https://github.com/nismod/jamaica-infrastructure">
        github.com/nismod/jamaica-infrastructure
      </ExtLink>
      .
    </p>
    <h2 id="navigate-around-the-j-srat">Navigate around the J-SRAT</h2>
    <a href="#table-of-contents">
      <small>Back to top</small>
    </a>
    <p>The top level of navigation for the tool is in the black top bar.</p>
    <p>“J-SRAT” brings you to the home page, which gives a brief introduction and summary of the analysis.</p>
    <p>
      Click across the links in the top navigation bar to see the “Exposure”, “Risk” and “Adaptation” stages of the
      flooding, cyclone and drought infrastructure risk assessment and adaptation analysis.
    </p>
    <p>“Nature-based Solutions” includes information about land-use and nature-based solutions.</p>
    <p>“Data” includes a summary of data used in the tool.</p>
    <p>
      <img src="/guide_media_media/image2.png" />
    </p>
    <p>Click on “Exposure” to show the map view. The main controls on this screen are used throughout.</p>
    <p>The left sidebar has various sections which control the data that is shown on the map.</p>
    <p>
      <img src="/guide_media/image3.png" />
    </p>
    <p>
      Click the search icon which is just to the right of the sidebar sections to search for places. This uses the
      OpenStreetMap Nominatim service and should find parishes, towns, some roads and some addresses by name. For
      example, search for “Norman Manley International”.
    </p>
    <p>
      <img src="/guide_media/image4.png" />
    </p>
    <p>Click the first result to zoom to the airport.</p>
    <p>
      <img src="/guide_media/image5.png" />
    </p>
    <p>Below the search box, there is a map layer control. Hover over the layers icon to show it.</p>
    <p>
      Switch the map background from the light grey “Map” background (designed by CARTO using OpenStreetMap data) to the
      “Satellite” imagery background (produced by EOX from Copernicus Sentinel-2 data).
    </p>
    <p>Check or uncheck the box to hide or “Show labels”</p>
    <p>
      <img src="/guide_media/image6.png" />
    </p>
    <p>
      In the top-right corner of the map, the plus and minus buttons control the map zoom. You can also scroll to zoom
      or double-click to zoom in and hold the shift key and double-click to zoom out.
    </p>
    <p>
      In the bottom-right corner of the map, there is a scale bar for reference and an “i” icon which shows or hides
      information about the background maps when clicked.
    </p>
    <h2 id="explore-the-j-srat-data-and-analysis">Explore the J-SRAT data and analysis</h2>
    <a href="#table-of-contents">
      <small>Back to top</small>
    </a>
    <h2 id="how-are-the-transport-energy-and-water-systems-represented">
      How are the transport, energy and water systems represented?
    </h2>
    <p>
      The “Exposure”, “Risk”, “Adaptation” and “Nature-based Solutions” tabs follow a consistent layout. The left
      sidebar controls the data that is shown on the map.
    </p>
    <p>Start on the “Exposure” tab.</p>
    <p>Click on the “Infrastructure” section to expand or collapse it.</p>
    <p>Use the “eye” icons to hide or show all of a section’s layers at once.</p>
    <p>
      Under “Infrastructure”, expand “Power” and “Transmission” and select “High Voltage” lines. Explore the other power
      system layers to see the elements of the network as they are included in the analysis.
    </p>
    <p>
      <img src="/guide_media/image7.png" />
    </p>
    <p>Use the top “Power” checkbox to deselect all power assets, and expand “Transport” to bring in roads layers.</p>
    <p>
      <img src="/guide_media/image8.png" />
    </p>
    <p>Similarly, under “Water”, bring in the water supply, irrigation and wastewater systems.</p>
    <p>
      <img src="/guide_media/image9.png" />
    </p>
    <h2 id="how-are-flooding-and-hurricanes-represented">How are flooding and hurricanes represented?</h2>
    <p>
      Open the “Hazards” section and select “River Flooding” to show potential flooding on the map. This is a return
      period map, showing the depth of flooding in any location across the island which is expected to be exceeded once
      in some number of years (once in 20 years, for example). Hover over shaded blue areas to see the depth of flooding
      in metres.
    </p>
    <p>
      <img src="/guide_media/image10.png" />
    </p>
    <p>
      In the sidebar, move the slider to show flooding for different return periods. A 500-year return period flood is
      much less likely and more intense, with deeper water levels and more area covered by the flood.
    </p>
    <p>
      <img src="/guide_media/image11.png" />
    </p>
    <p>While looking at the hazards, we can overlay infrastructure networks to see where they might be affected.</p>
    <p>
      <img src="/guide_media/image12.png" />
    </p>
    <p>
      Click on a road, for example, to see details of the asset damage calculated from the length of road exposed to
      different depths of flooding at each return period.
    </p>
    <p>
      <img src="/guide_media/image13.png" />
    </p>
    <h2 id="how-are-people-and-economic-activities-represented">How are people and economic activities represented?</h2>
    <p>Population is mapped to administrative boundaries, and economic activity is assigned to buildings.</p>
    <p>
      In the left sidebar, open the “Regions” tab to hide or show boundaries. “Parishes” shows the top-level regions.
    </p>
    <p>
      <img src="/guide_media/image14.png" />
    </p>
    <p>Select “Enumeration Districts” to see the small areas.</p>
    <p>
      <img src="/guide_media/image15.png" />
    </p>
    <p>
      Change the “Layer Style” to “Population” to show population density on the map. Hover over areas to see population
      counts, and click on an area for a small detail sidebar to appear on the right of the screen.
    </p>
    <p>
      <img src="/guide_media/image16.png" />
    </p>
    <p>
      Expand the “Buildings” section in the left sidebar and check that the eye icon is toggled on, then zoom in to see
      buildings. The buildings are not shown at all until quite high zoom levels.
    </p>
    <p>
      <img src="/guide_media/image17.png" />
    </p>
    <p>
      Click on a building to see details, including the total assigned GDP and estimated rehabilitation cost. Risk and
      direct damages from flooding and cyclones have been assessed for buildings as for infrastructure assets, though no
      indirect effects are estimated beyond the disruption to the activity in the building itself.
    </p>
    <p>
      <img src="/guide_media/image18.png" />
    </p>
    <p>Toggle the “Hazards” section visibility to show flood maps under buildings.</p>
    <p>
      <img src="/guide_media/image19.png" />
    </p>
    <h2 id="how-is-climate-risk-from-flooding-or-cyclones-represented">
      How is climate risk from flooding or cyclones represented?
    </h2>
    <p>
      On the “Risk” page, in the “Infrastructure” sidebar section, the “Layer style” defaults to show damages. Select
      the “Power”, “Transmission”, “Substations” layer to show expected direct damages or economic losses as evaluated
      for all of the substations across the island.
    </p>
    <p>
      <img src="/guide_media/image20.png" />
    </p>
    <p>
      From the controls in the sidebar, explore the contribution of individual hazards, and select the epoch (year) and
      RCP (climate scenario) to see how risks change as the hazards change. At this point in the analysis there is no
      accounting for economic growth or discounting over time: the comparison is holding everything constant except for
      the hazard maps.
    </p>
    <h2 id="how-is-drought-risk-represented">How is drought risk represented?</h2>
    <p>
      On the “Adaptation” page, the second sidebar section contains controls to show drought risk and drought-related
      adaptation options. Select for example “Population at risk” to see areas coloured by risk and "Population
      protected" to see points representing approximate locations of options to reduce the impacts of drought.
    </p>
    <p>
      <img src="/guide_media/image21.png" />
    </p>
    <h2 id="how-are-adaptation-options-presented">How are adaptation options presented?</h2>
    <p>
      On the “Adaptation” page, the left sidebar includes an “Adaption Options” layer style. Select a sector, sub sector
      and asset type to display the generic adaptation options that have been evaluated. For example, for “Transport”,
      “Road”, “Class A”, you can look at “Elevate the roads” or a combined upgrade option.
    </p>
    <p>
      <img src="/guide_media/image22.png" />
    </p>
    <p>
      From here, the assets are shown on the map and ranked in the table on the right, according to the choice of
      “Displayed variable”, which includes net present value of avoided direct damages or economic losses, total
      adaptation cost, or the cost-benefit ratio.
    </p>
    <p>
      Hover over the rows of the table on the right to highlight asset locations. Click the magnifying glass icon to
      zoom in. Click the “zoom out” magnifying glass icon at the top of the table to get back to the full island view.
      Click on a row of the table to see the main asset attributes.
    </p>
    <p>
      Once you have identified an asset of interest, click on it on the map to highlight it in blue, the switch the
      layer style to “Damages” – or switch the page to “Risk” – to show the full asset details sidebar. Scroll down for
      the full details on Adaptation Options for this particular asset.
    </p>
    <p>
      <img src="/guide_media/image23.png" />
    </p>
    <h2 id="how-are-nature-based-solutions-represented">How are nature-based solutions represented?</h2>
    <p>The “Nature-based solutions” tab enables analysis and exploration of various aspects of the land and sea.</p>
    <p>Expand the “Terrestrial” section to show a long list of Land Use/Land Cover classes.</p>
    <p>
      <img src="/guide_media/image24.png" />
    </p>
    <p>
      Change the maximum or minimum values for elevation or slope to constrain the areas displayed – for example, look
      for high-slope areas with land cover that could be improved for slope stability and runoff reduction.
    </p>
    <p>Apply various other constraints to see areas which are protected, or within 100m of a stream or forest.</p>
    <p>
      Expand the “Marine” section to show coral, mangrove and seagrass, with 500m buffer zones around each habitat area.
    </p>
    <p>
      <img src="/guide_media/image25.png" />
    </p>
    <h2 id="trace-the-climate-risk-and-adaptation-analysis-for-a-single-asset">
      Trace the climate risk and adaptation analysis for a single asset
    </h2>
    <a href="#table-of-contents">
      <small>Back to top</small>
    </a>
    <p>
      This use case revisits many of the previous stages, looking at a single infrastructure asset to understand the
      process of climate risk and adaptation analysis. It assumes that there is a particular asset that you are
      interested in and you want to assess it for risk and potential adaptation.
    </p>
    <p>This example will use the substation with id “node_17”, near Porus in the Parish of Manchester.</p>
    <p>Start on the “Exposure” tab.</p>
    <p>
      Under “Infrastructure” in the left sidebar, select “Power &gt; Transmission &gt; Substations” to show all
      substations.
    </p>
    <p>In the map search box, type “Porus” and click on “Porus, Manchester” to zoom to the town and substation.</p>
    <p>Under “Hazards” in the left sidebar, select “River flooding” to add a flood outline to the map.</p>
    <p>Hover over the asset to see its ID and depth of flooding.</p>
    <p>
      <img src="/guide_media/image26.png" />
    </p>
    <p>
      The substation intersects with different fluvial (river) flood outlines under baseline and future climate
      scenarios. From the depth of flooding, the analysis has estimated the direct damages to the substation itself, and
      indirect economic losses from buildings and businesses experiencing loss of power.
    </p>
    <p>
      Click on the substation to show the calculated results in the right-hand details sidebar. Under the “Risk” section
      (not shown), the chart and table show the expected annual damages or losses.
    </p>
    <p>
      Under the “Return Period Damages” section (shown below), the chart and table display in more detail the damages or
      losses which would result from a flood of a particular return period (or probability).
    </p>
    <p>
      <img src="/guide_media/image27.png" />
    </p>
    <p>
      In the “Return Period Damages” section, change the “Epoch” dropdown from 2010 to 2080 to see the change in damages
      under various possible climate change scenarios (i.e. Representative Concentration Pathways [RCP] 2.6, 4.5 and
      8.5).
    </p>
    <p>
      <img src="/guide_media/image28.png" />
    </p>
    <p>
      In the top menu, click on the “Risk” tab. This may change the background flood map but should leave the map
      location and asset selection unchanged.
    </p>
    <p>Scroll back up to the “Risk” section in the right-hand sidebar.</p>
    <p>
      <img src="/guide_media/image29.png" />
    </p>
    <p>
      One important thing to note is that (so far) only the flood hazard part of the equation has been changing in the
      calculations above. The “2080” risk values and return period damages have not (yet) considered any changes in
      population or the economy. These changes are included in the next step.
    </p>
    <p>
      Once all the climate risks have been estimated, the analysis considers a chosen adaptation option for this
      substation, which is to build a protective wall that prevents flood waters from inundating the asset. We test
      different heights of flood protection walls and estimate the Net Present Value (NPV) cost, NPV benefits (in terms
      of avoided risks), and Benefit-Cost Ratio (BCR) values for each case.
    </p>
    <p>
      Net Present Value is calculated over the asset lifetime, and both economic growth scenarios and discount rates go
      into the calculation. Population change is considered in some sectors – particularly in water, when calculating
      the risk of drought. The adaptation investment and maintenance costs, or the direct or indirect risk values are
      calculated for each year into the future, discounted and summed to give a single Net Present Value number. For
      more detail on these calculations, see section 2.2.2 of the Technical Report.
    </p>
    <p>Scroll down to the “Adaptation Options” section.</p>
    <p>
      <img src="/guide_media/image30.png" />
    </p>
    <p>
      Here we see that building a protective wall that is 1-meter high would incur an NPV cost of J$ 110 million over
      time and result in avoiding mean NPV risks of J$ 807 (somewhere between 761 and 863) million under RCP 2.6
      flooding scenario, which results in BCR of 7.32 (somewhere between 6.91 and 7.83), which is much greater than the
      break-even value of 1. For this asset the chosen adaptation option presents a robust case for investment: the
      benefits are much greater than the costs.
    </p>
    <h2 id="prioritise-and-evaluate-adaptation-interventions">Prioritise and evaluate adaptation interventions</h2>
    <a href="#table-of-contents">
      <small>Back to top</small>
    </a>
    <p>
      This use case covers reviewing and prioritising adaptation options for an infrastructure sector, then going into
      detail on a specific asset to evaluate the direct and indirect risks, and costs and benefits of adaptation. It
      assumes that you start from a system- or sector- wide perspective and want to go through a process of screening
      assets for risk and potential adaptation.
    </p>
    <p>Start from the “Adaptation” tab, chosen from the top menu bar.</p>
    <p>In the left sidebar, make sure “Adaptation Options” is the layer style.</p>
    <p>
      Select the sector as “Power”, subsector as “Transmission” and asset type as “Substation” using the dropdown menus
      (the sector check boxes on top left are not functional on this tab).
    </p>
    <p>
      Select details of the hazard (“Flooding” is the only one available for this asset type) against which the
      adaptation option will protect. Here you can also change the climate scenario, adaptation option type, and (for
      some option types) protection level.
    </p>
    <p>Select “Cost-Benefit Ratio” as the displayed variable.</p>
    <p>
      The map then shows assets coloured by Cost-Benefit Ratio for prioritisation. In the right sidebar, the table shows
      assets, sorted in descending order with the most cost-beneficial at the top.
    </p>
    <p>
      <img src="/guide_media/image31.png" />
    </p>
    <p>Hover over a row to indicate the asset on the map, drawing a dashed bright blue line around its location.</p>
    <p>
      Click on a row for a few more details about the asset. At the right-hand end of the row, click on the zoom-in
      magnifying glass icon (🔍 ) with a plus sign to zoom to the asset location (
      <em>you will need to move the cursor back onto the row to see this, if it has moved elsewhere</em>).
    </p>
    <p>
      To zoom out again to the whole island, click on the magnifying glass icon (🔍 ) with a minus sign, which is at the
      top-right corner of the sidebar.
    </p>
    <p>
      <img src="/guide_media/image32.png" />
    </p>
    <p>
      Once you have identified a candidate for prioritisation, zoom to its location so that it is clearly in view. You
      can then switch back to the “Exposure” or “Risk” tabs to bring up more details about the asset risk. In this
      example, we have prioritised the same substation with id “node_17”, near Porus.
    </p>
    <p>
      Click on the “Risk” tab in the top menu bar, then click on the asset to highlight it and show the right-hand
      details sidebar. Here, as we have seen previously, there are details of asset risk under different hazards and
      climate scenarios, both direct damages and indirect economic losses.
    </p>
    <p>
      <img src="/guide_media/image33.png" />
    </p>
    <p>
      Scroll down to the “Adaptation Options” section to find the evaluated Net Present Value costs and benefits of
      different adaptation interventions under different hazards and climate scenarios.
    </p>
    <p>
      <img src="/guide_media/image34.png" />
    </p>
    <p>
      From any of the asset details sections (e.g. Adaptation Options as shown above), click on the download icon (⤓) to
      save a CSV of all the table values, including all hazards/epochs.
    </p>

    <h2 id="about-this-guide">About this guide</h2>
    <a href="#table-of-contents">
      <small>Back to top</small>
    </a>
    <p>This document may be cited as follows:</p>
    <blockquote>
      Russell, Tom, Ziarkowski, Maciej, Pant, Raghav, Becher, Olivia, Haggis Robyn, Majid, Aman, Verschuur, Jasper,
      Grant, Ardith, Mills, Anaitee, Williams, Edson, Fowler, Tim and Hall, Jim W. (2022).{' '}
      <em>User Guide for the Jamaica Systemic Risk Assessment Tool (J-SRAT)</em>. Environmental Change Institute,
      University of Oxford, UK.
    </blockquote>
    <p>Report © University of Oxford, 2022</p>
    <p>
      This report was produced as part of the project “A geospatial analysis platform for infrastructure risk assessment
      and resilient investment prioritisation in Jamaica” which was funded by the UK Foreign Commonwealth and
      Development Office to inform the (FDCO) as part of the Coalition for Climate Resilient Investments. The views
      expressed and recommendations set out in this report are the authors’ own and do not necessarily reflect the
      position of the FCDO, CCRI or any stakeholder in Jamaica.
    </p>
    <p>
      The materials have been prepared by Oxford University. Whilst every care has been taken by Oxford University to
      ensure the accuracy and completeness of the reports and maps, the reader must recognise that errors are possible
      through no fault of Oxford University and as such the parties give no express or implied representations or
      warranty as to:
    </p>
    <p>
      (i) the quality or fitness for any particular purpose of the report or maps supplied or of any design,
      workmanship, materials or parts used in connection therewith or correspondence with regard to any description or
      sample; or
    </p>
    <p>
      (ii) the accuracy, sufficiency or completeness of the reports or maps provided. In particular, there are hereby
      expressly excluded all conditions, warranties and other terms which might otherwise be implied (whether by common
      law, by statute or otherwise).
    </p>
    <p>
      Oxford University, its employees, servants and agents shall accept no liability for any damage caused directly or
      indirectly by the use of any information contained herein and without prejudice to the generality of the
      foregoing, by any inaccuracies, defects or omissions.
    </p>
  </article>
);
