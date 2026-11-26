# Hollinger

<div align = "center">
<img src = "https://i.postimg.cc/3ND2B1zL/crd-coat-of-arms-logo.png" height = "64">
<br><br>

https://github.com/user-attachments/assets/9a571209-1e1f-4ce8-a58c-b1dc5a552737

<div align = "center">
  Figure 1. Global land surface change and deglaciation from 10000BC-2025AD at annual resolution.
</div>
<br>

[![Join our community (Element Matrix)!](https://img.shields.io/badge/chat-on%20matrix-51bb9c?style=for-the-badge)](https://matrix.to/#/#confoederatio:matrix.confoederatio.org) [![Join our community!](https://img.shields.io/discord/548994743925997570?label=Discord&style=for-the-badge)](https://discord.gg/89kQY2KFQz) ![](https://img.shields.io/github/languages/code-size/ConfoederatioVF/Velkscala?style=for-the-badge)

</div>

## Abstract.

- E-mail: [vf@confoederatio.org](mailto:vf@confoederatio.org)
- Collation: [Github](https://github.com/ConfoederatioVF/Collation) | [Projects](https://projects.confoederatio.org/c/Collation)
- Dataview: https://dataview.confoederatio.org/
- <ins>Download:</ins>
- - Stadestér (Cities, part of Velkscala) | [Github](https://github.com/ConfoederatioVF/Stadester) | [Projects (1.1)](https://projects.confoederatio.org/c/Stadester)
- - **Velkscala (Demographics)** | [Github](https://github.com/ConfoederatioVF/Velkscala) | [Datavault (1.0)](https://datavault.confoederatio.org/s/44Pfgt7ya8wE77b) - Includes <ins>Hollinger</ins>.
- - Eoscala (Economics) | [Github](https://github.com/ConfoederatioVF/Eoscala) | [Datavault (1.5)](https://datavault.confoederatio.org/s/2Kytdej7ABpw4rr)
- - Atlas (Polities) | [Github](https://github.com/ConfoederatioVF/Atlas) | [Datavault (0.7)](https://datavault.confoederatio.org/s/rPntqXwCxzXg4Bk)

**Hollinger** is a land cover and land use change dataset from 10000BC-2025AD, meant as a drop-in replacement for HYDE3.5. Because HYDE does not factor in coastline change and small island nations or marginal coastlines, the commissioning of a separate land use estimate dataset was necessary for priors to be as accurate as possible.

As a sub-dataset of [Velkscala](https://github.com/ConfoederatioVF/Velkscala), Hollinger is a part of the **Naissance World Model**, a digital model of human society from 10000BC-2025AD.

## File Structure.

All rasters are rendered at 4320x2160 (WGS84 Equirectangular), with 100-year timesteps from 10000BC-1700AD (where 0 is 1AD), 10-year timesteps from 1700AD-1950AD, and annual timesteps since to 2025AD.

**Metadata (float32):**
- `landarea.nc`: Land/Glacier/Sea surface area (km^2). Land is positive, glaciers are negative (abs. value), and sea is zero.

**Biomes & Climate (int16):**

All biome and climate NetCDF files are encoded as integer classification schemes.

- `anthromes.nc`: Replacement Anthromes12k classifications.
  - 0: Water
  - 11: Urban; 12: Mixed settlements
  - 21: Rice villages; 22: Irrigated villages, 23: Rainfed villages, 24: Pastoral villages
  - 31: Residential irrigated croplands, 32: Residential rainfed croplands, 33: Populated croplands, 34: Remote croplands
  - 41: Residential rangelands, 42: Populated rangelands, 43: Remote rangelands
  - 51: Residential woodlands, 52: Populated woodlands, 53: Remote woodlands, 54: Inhabited drylands
  - 61: Wild woodlands, 62: Wild drylands, 63: Ice
- `holdridge_lifezones.nc`: Holdridge lifezones (Biomes).
  - 0: Water
  - 1: Polar desert
  - 2: Subpolar dry tundra, 3: Subpolar moist tundra, 4: Subpolar wet tundra, 5: Subpolar rain tundra
  - 6: Boreal desert, 7: Boreal dry scrub, 8: Boreal moist forest, 9: Boreal wet forest, 10: Boreal rain forest
  - 11: Cool temperate desert, 12: Cool temperate desert scrub, 13: Cool temperate steppe, 14: Cool temperate moist forest, 15: Cool temperate wet forest, 16: Cool temperate rain forest
  - 17: Warm temperate desert, 18: Warm temperate desert scrub, 19: Warm temperate thorn steppe, 20: Warm temperate dry forest, 21: Warm temperate moist forest, 22: Warm temperate wet forest, 23: Warm temperate rain forest
  - 24: Subtropical desert, 25: Subtropical desert scrub, 26: Subtropical thorn woodland, 27: Subtropical dry forest, 28: Subtropical moist forest, 29: Subtropical wet forest, 30: Subtropical rain forest
  - 31: Tropical desert, 32: Tropical desert scrub, 33: Tropical thorn woodland, 34: Tropical very dry forest, 35: Tropical dry forest, 36: Tropical moist forest, 37: Tropical wet forest, 38: Tropical rain forest
- `koeppen_geiger.nc`: Köppen-Geiger (Climate).
  - 0: Water
  - 1: Af (Tropical, rainforest), 2: Am (Tropical, monsoon), 3: Aw (Tropical, savannah)
  - 4: BWh (Arid, desert, hot), 5: BWk (Arid, desert, cold), 6: BSh (Arid, steppe, hot), 7: BSk (Arid, steppe, cold)
  - 8: Csa (Temperate, dry summer, hot summer), 9: Csb (Temperate, dry summer, warm summer), 10: Csc (Temperate, dry summer, cold summer)
  - 11: Cwa (Temperate, dry winter, hot summer), 12: Cwb (Temperate, dry winter, warm summer), 13: Cwc (Temperate, dry winter, cold summer)
  - 14: Cfa (Temperate, no dry season, hot summer), 15: Cfb (Temperate, no dry season, warm summer), 16: Cfc (Temperate, no dry season, cold summer)
  - 17: Dsa (Cold, dry summer, hot summer), 18: Dsb (Cold, dry summer, warm summer), 19: Dsc (Cold, dry summer, cold summer), 20: Dsd (Cold, dry summer, very cold winter)
  - 21: Dwa (Cold, dry winter, hot summer), 22: Dwb (Cold, dry winter, warm summer), 23: Dwc (Cold, dry winter, cold summer), 24: Dwd (Cold, dry winter, very cold winter)
  - 25: Dfa (Cold, no dry season, hot summer), 26: Dfb (Cold, no dry season, warm summer), 27: Dfc (Cold, no dry season, cold summer), 28: Dfd (Cold, no dry season, very cold winter)
  - 29: ET (Polar, tundra), 30: EF (Polar, frost)
- `trewartha.nc`: Trewartha (Climate).
  - 0: Water
  - 1: Ar (Tropical wet), 2: Am (Tropical monsoon), 3: Aw (Tropical wet-dry, dry winter), 4: As (Tropical wet-dry, dry summer)
  - 5: BWh (Subtropical/tropical desert), 6: BWk (Temperate/boreal desert), 7: BSh (Subtropical/tropical steppe), 8: BSk (Temperate/boreal steppe)
  - 9: Cs (Subtropical dry summer), 10: Cw (Subtropical dry winter), 11: Cr (Subtropical humid)
  - 12: Do (Temperate oceanic), 13: Dc (Temperate continental)
  - 14: Eo (Boreal oceanic), 15: Ec (Boreal continental)
  - 16: Ft (Polar tundra), 17: Fi (Polar ice cap)

## Sources/Attribution.

<details>
  <summary>Anthromes (Expand to view sources).</summary>

1. Ellis, Erle C., et al. “Anthropogenic Biomes: 10,000 BCE to 2015 CE.” *Land*, vol. 9, no. 5, 2020, 129, <https://doi.org/10.3390/land9050129>.
2. Ellis, Erle C., et al. “People Have Shaped Most of Terrestrial Nature for at Least 12,000 Years.” *Proceedings of the National Academy of Sciences*, vol. 118, no. 17, 2021, e2023483118, <https://doi.org/10.1073/pnas.2023483118>.
3. Gauthier, Nicolas, et al. *Anthromes 12K DGG (V1) Full Dataset*. Harvard Dataverse, 2021, <https://doi.org/10.7910/DVN/E3H3AK>.
4. Ramankutty, Navin, and Jonathan A. Foley. “Estimating Historical Changes in Global Land Cover: Croplands from 1700 to 1992.” *Global Biogeochemical Cycles*, vol. 13, no. 4, 1999, pp. 997-1027, <https://doi.org/10.1029/1999GB900046>.
5. Klein Goldewijk, Kees. *HYDE Data for the Global Carbon Budget (2025)*. Utrecht University - Yoda, 2025, <https://doi.org/10.24416/UU01-F45D44>.
6. Klein Goldewijk, Kees, et al. “Anthropogenic Land Use Estimates for the Holocene - HYDE 3.2.” *Earth System Science Data*, vol. 9, no. 2, 2017, pp. 927-953, <https://doi.org/10.5194/essd-9-927-2017>.
7. Schiavina, Marcello, et al. *GHS-POP R2023A - GHS Population Grid Multitemporal (1975-2030)*. European Commission, Joint Research Centre, 2023, <https://doi.org/10.2905/2FF68A52-5B5B-4A22-8F40-C41DA8332CFE>.
8. Melchiorri, Michele, et al. “The Multi-Temporal and Multi-Dimensional Global Urban Centre Database to Delineate and Analyse World Cities.” *Scientific Data*, vol. 11, no. 1, 2024, 82, <https://doi.org/10.1038/s41597-023-02691-1>.
9. Marí Rivero, Inés, et al. *GHS-UCDB R2024A - GHS Urban Centre Database 2025*. European Commission, Joint Research Centre, 2024, <https://doi.org/10.2905/1A338BE6-7EAF-480C-9664-3A8ADE88CBCD>.
10. Gapminder Foundation. *Gapminder Population Dataset (Version 6)*. Gapminder, 2022, <https://www.gapminder.org/data/documentation/gd003/>.
  
</details>

<details>
  <summary>Holdridge Lifezones (Expand to view sources).</summary>

1. Holdridge, Leslie R. “Determination of World Plant Formations from Simple Climatic Data.” *Science*, vol. 105, no. 2727, 1947, pp. 367-368, <https://doi.org/10.1126/science.105.2727.367>.
2. Holdridge, Leslie R. *Life Zone Ecology*. Rev. ed., Tropical Science Center, 1967, <https://search.worldcat.org/search?q=Life+Zone+Ecology+Holdridge+1967>.
3. Leemans, Rik. *Possible Changes in Natural Vegetation Patterns Due to Global Warming*. IIASA Working Paper WP-90-08, International Institute for Applied Systems Analysis, 1990, <https://pure.iiasa.ac.at/id/eprint/3443/>.
4. Lugo, Ariel E., et al. “The Holdridge Life Zones of the Conterminous United States in Relation to Ecosystem Mapping.” *Journal of Biogeography*, vol. 26, no. 5, 1999, pp. 1025-1038, <https://doi.org/10.1046/j.1365-2699.1999.00329.x>.
5. Liu, Zhengyu, et al. “Transient Simulation of Last Deglaciation with a New Mechanism for Bølling-Allerød Warming.” *Science*, vol. 325, no. 5938, 2009, pp. 310-314, <https://doi.org/10.1126/science.1171041>.
6. Otto-Bliesner, Bette, and Nan Rosenbloom. *Simulation of the Transient Climate of the Last 21,000 Years (TraCE-21ka)*. NSF National Center for Atmospheric Research, 2021, <https://doi.org/10.5065/CXB5-TV56>.
7. Fick, Stephen E., and Robert J. Hijmans. “WorldClim 2: New 1-km Spatial Resolution Climate Surfaces for Global Land Areas.” *International Journal of Climatology*, vol. 37, no. 12, 2017, pp. 4302-4315, <https://doi.org/10.1002/joc.5086>.
8. WorldClim. *WorldClim Version 2.1 Historical Climate Data (1970-2000, 5 Arc-Minute)*. University of California, Davis, 2020, <https://www.worldclim.org/data/worldclim21.html>.
9. Harris, Ian, et al. “Version 4 of the CRU TS Monthly High-Resolution Gridded Multivariate Climate Dataset.” *Scientific Data*, vol. 7, no. 1, 2020, 109, <https://doi.org/10.1038/s41597-020-0453-3>.
10. Climatic Research Unit, University of East Anglia. *CRU TS v4.10 Gridded Time-Series Dataset (1901-2025)*. University of East Anglia, 2026, <https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.10/>.
11. Tarasov, Lev, et al. “A Data-Calibrated Distribution of Deglacial Chronologies for the North American Ice Complex from Glaciological Modeling.” *Earth and Planetary Science Letters*, vol. 315-316, 2012, pp. 30-40, <https://doi.org/10.1016/j.epsl.2011.09.010>.
12. Briggs, Robert D., et al. “A Data-Constrained Large Ensemble Analysis of Antarctic Evolution since the Eemian.” *Quaternary Science Reviews*, vol. 103, 2014, pp. 91-115, <https://doi.org/10.1016/j.quascirev.2014.09.003>.
13. Ivanovic, Ruza F., et al. “Transient Climate Simulations of the Deglaciation 21-9 Thousand Years before Present (Version 1) - PMIP4 Core Experiment Design and Boundary Conditions.” *Geoscientific Model Development*, vol. 9, no. 7, 2016, pp. 2563-2587, <https://doi.org/10.5194/gmd-9-2563-2016>.
  
</details>

<details>
  <summary>Köppen-Geiger (Expand to view sources).</summary>

1. Peel, Murray C., et al. “Updated World Map of the Köppen-Geiger Climate Classification.” *Hydrology and Earth System Sciences*, vol. 11, no. 5, 2007, pp. 1633-1644, <https://doi.org/10.5194/hess-11-1633-2007>.
2. Beck, Hylke E., et al. “Present and Future Köppen-Geiger Climate Classification Maps at 1-km Resolution.” *Scientific Data*, vol. 5, no. 1, 2018, 180214, <https://doi.org/10.1038/sdata.2018.214>.
3. Beck, Hylke E., et al. “High-Resolution (1 km) Köppen-Geiger Maps for 1901-2099 Based on Constrained CMIP6 Projections.” *Scientific Data*, vol. 10, no. 1, 2023, 724, <https://doi.org/10.1038/s41597-023-02549-6>.
4. Liu, Zhengyu, et al. “Transient Simulation of Last Deglaciation with a New Mechanism for Bølling-Allerød Warming.” *Science*, vol. 325, no. 5938, 2009, pp. 310-314, <https://doi.org/10.1126/science.1171041>.
5. Otto-Bliesner, Bette, and Nan Rosenbloom. *Simulation of the Transient Climate of the Last 21,000 Years (TraCE-21ka)*. NSF National Center for Atmospheric Research, 2021, <https://doi.org/10.5065/CXB5-TV56>.
6. Fick, Stephen E., and Robert J. Hijmans. “WorldClim 2: New 1-km Spatial Resolution Climate Surfaces for Global Land Areas.” *International Journal of Climatology*, vol. 37, no. 12, 2017, pp. 4302-4315, <https://doi.org/10.1002/joc.5086>.
7. WorldClim. *WorldClim Version 2.1 Historical Climate Data (1970-2000, 5 Arc-Minute)*. University of California, Davis, 2020, <https://www.worldclim.org/data/worldclim21.html>.
8. Harris, Ian, et al. “Version 4 of the CRU TS Monthly High-Resolution Gridded Multivariate Climate Dataset.” *Scientific Data*, vol. 7, no. 1, 2020, 109, <https://doi.org/10.1038/s41597-020-0453-3>.
9. Climatic Research Unit, University of East Anglia. *CRU TS v4.10 Gridded Time-Series Dataset (1901-2025)*. University of East Anglia, 2026, <https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.10/>.
10. Tarasov, Lev, et al. “A Data-Calibrated Distribution of Deglacial Chronologies for the North American Ice Complex from Glaciological Modeling.” *Earth and Planetary Science Letters*, vol. 315-316, 2012, pp. 30-40, <https://doi.org/10.1016/j.epsl.2011.09.010>.
11. Briggs, Robert D., et al. “A Data-Constrained Large Ensemble Analysis of Antarctic Evolution since the Eemian.” *Quaternary Science Reviews*, vol. 103, 2014, pp. 91-115, <https://doi.org/10.1016/j.quascirev.2014.09.003>.
12. Ivanovic, Ruza F., et al. “Transient Climate Simulations of the Deglaciation 21-9 Thousand Years before Present (Version 1) - PMIP4 Core Experiment Design and Boundary Conditions.” *Geoscientific Model Development*, vol. 9, no. 7, 2016, pp. 2563-2587, <https://doi.org/10.5194/gmd-9-2563-2016>.
  
</details>

<details>
  <summary>Land area (Expand to view sources).</summary>

1. OpenStreetMap Contributors. *Land Polygons (Split, WGS84)*. FOSSGIS e.V., 2025, <https://osmdata.openstreetmap.de/data/land-polygons.html>.
2. NOAA National Centers for Environmental Information. *ETOPO 2022 15 Arc-Second Global Relief Model*. NOAA NCEI, 2022, <https://doi.org/10.25921/fd45-gt74>.
3. Tarasov, Lev, et al. “A Data-Calibrated Distribution of Deglacial Chronologies for the North American Ice Complex from Glaciological Modeling.” *Earth and Planetary Science Letters*, vol. 315-316, 2012, pp. 30-40, <https://doi.org/10.1016/j.epsl.2011.09.010>.
4. Briggs, Robert D., et al. “A Data-Constrained Large Ensemble Analysis of Antarctic Evolution since the Eemian.” *Quaternary Science Reviews*, vol. 103, 2014, pp. 91-115, <https://doi.org/10.1016/j.quascirev.2014.09.003>.
5. Ivanovic, Ruza F., et al. “Transient Climate Simulations of the Deglaciation 21-9 Thousand Years before Present (Version 1) - PMIP4 Core Experiment Design and Boundary Conditions.” *Geoscientific Model Development*, vol. 9, no. 7, 2016, pp. 2563-2587, <https://doi.org/10.5194/gmd-9-2563-2016>.
6. Messager, Mathis Loïc, et al. “Estimating the Volume and Age of Water Stored in Global Lakes Using a Geo-Statistical Approach.” *Nature Communications*, vol. 7, 2016, 13603, <https://doi.org/10.1038/ncomms13603>.
7. Lehner, Bernhard, et al. “The Global Dam Watch Database of River Barrier and Reservoir Information for Large-Scale Applications.” *Scientific Data*, vol. 11, 2024, 1069, <https://doi.org/10.1038/s41597-024-03752-9>.
8. Lehner, Bernhard, et al. *Global Dam Watch (GDW) Database v1.0*. Figshare, 2024, <https://doi.org/10.6084/m9.figshare.25988293>.
9. Natural Earth. *1:10m Physical Vectors: Glaciated Areas; Antarctic Ice Shelves*. North American Cartographic Information Society, <https://www.naturalearthdata.com/downloads/10m-physical-vectors/>.
10. Lambeck, Kurt, et al. “Sea Level and Global Ice Volumes from the Last Glacial Maximum to the Holocene.” *Proceedings of the National Academy of Sciences*, vol. 111, no. 43, 2014, pp. 15296-15303, <https://doi.org/10.1073/pnas.1411762111>.
11. Peltier, W. Richard, et al. “Space Geodesy Constrains Ice Age Terminal Deglaciation: The Global ICE-6G_C (VM5a) Model.” *Journal of Geophysical Research: Solid Earth*, vol. 120, no. 1, 2015, pp. 450-487, <https://doi.org/10.1002/2014JB011176>.
12. Norðdahl, Hreggviður, and Hlynur G. Pétursson. “Relative Sea-Level Changes in Iceland: New Aspects of the Weichselian Deglaciation of Iceland.” *Developments in Quaternary Sciences*, vol. 5, Elsevier, 2005, pp. 25-78, <https://doi.org/10.1016/S1571-0866(05)80005-3>.
13. Patton, Henry, et al. “The Configuration, Sensitivity and Rapid Retreat of the Late Weichselian Icelandic Ice Sheet.” *Earth-Science Reviews*, vol. 166, 2017, pp. 223-245, <https://doi.org/10.1016/j.earscirev.2017.02.001>.
14. Davies, Bethan J., et al. “The Evolution of the Patagonian Ice Sheet from 35 ka to the Present Day (PATICE).” *Earth-Science Reviews*, vol. 204, 2020, 103152, <https://doi.org/10.1016/j.earscirev.2020.103152>.
15. Menounos, Brian, et al. “Cordilleran Ice Sheet Mass Loss Preceded Climate Reversals near the Pleistocene Termination.” *Science*, vol. 358, no. 6364, 2017, pp. 781-784, <https://doi.org/10.1126/science.aan3001>.
16. Dalton, April S., et al. “An Updated Radiocarbon-Based Ice Margin Chronology for the Last Deglaciation of the North American Ice Sheet Complex.” *Quaternary Science Reviews*, vol. 234, 2020, 106223, <https://doi.org/10.1016/j.quascirev.2020.106223>.
17. Dyke, Arthur S. “An Outline of North American Deglaciation with Emphasis on Central and Northern Canada.” *Developments in Quaternary Sciences*, vol. 2, Elsevier, 2004, pp. 373-424, <https://doi.org/10.1016/S1571-0866(04)80209-4>.
18. Barrell, David J. A. “Quaternary Glaciers of New Zealand.” *Developments in Quaternary Sciences*, vol. 15, Elsevier, 2011, pp. 1047-1064, <https://doi.org/10.1016/B978-0-444-53447-7.00075-1>.
19. Bentley, Michael J., et al. “A Community-Based Geological Reconstruction of Antarctic Ice Sheet Deglaciation since the Last Glacial Maximum.” *Quaternary Science Reviews*, vol. 100, 2014, pp. 1-9, <https://doi.org/10.1016/j.quascirev.2014.06.025>.
20. Vos, Peter, et al. *Atlas of the Holocene Netherlands: Landscape and Habitation since the Last Ice Age*. Amsterdam University Press, 2020, <https://doi.org/10.5117/9789463724432>.
21. Rijksdienst voor het Cultureel Erfgoed. *Paleogeografische Kaarten (WFS)*. RCE, <https://services.rce.geovoorziening.nl/paleogeografischekaarten/wfs>.
22. Hijma, Marc P., and Kim M. Cohen. “Timing and Magnitude of the Sea-Level Jump Preluding the 8200 yr Event.” *Geology*, vol. 38, no. 3, 2010, pp. 275-278, <https://doi.org/10.1130/G30439.1>.
23. Teller, James T., and David W. Leverington. “Glacial Lake Agassiz: A 5000 yr History of Change and Its Relationship to the δ18O Record of Greenland.” *Geological Society of America Bulletin*, vol. 116, no. 5, 2004, pp. 729-742, <https://doi.org/10.1130/B25316.1>.
24. Leverington, David W., et al. “Changes in the Bathymetry and Volume of Glacial Lake Agassiz between 9200 and 7700 14C yr B.P.” *Quaternary Research*, vol. 57, no. 2, 2002, pp. 244-252, <https://doi.org/10.1006/qres.2001.2311>.
25. Barber, Donald C., et al. “Forcing of the Cold Event of 8,200 Years Ago by Catastrophic Drainage of Laurentide Lakes.” *Nature*, vol. 400, no. 6742, 1999, pp. 344-348, <https://doi.org/10.1038/22504>.
26. Krijgsman, Wout, et al. “Quaternary Time Scales for the Pontocaspian Domain: Interbasinal Connectivity and Faunal Evolution.” *Earth-Science Reviews*, vol. 188, 2019, pp. 1-40, <https://doi.org/10.1016/j.earscirev.2018.10.013>.
27. Kuzucuoğlu, Catherine, and Suzanne A. G. Leroy. “Geographic and Geomorphologic Context of Caspian Sea Level Fluctuations over the Late Pleistocene and Holocene.” *Quaternaire*, vol. 34, no. 2, 2023, pp. 71-92, <https://doi.org/10.4000/quaternaire.17991>.
28. Svitoch, Alexander. “Holocene Transgressions of the Caspian Sea.” *Quaternary International*, vol. 279-280, 2012, pp. 478-479, <https://doi.org/10.1016/j.quaint.2012.08.1616>.
29. Yanko-Hombach, Valentina, et al., editors. *The Black Sea Flood Question: Changes in Coastline, Climate, and Human Settlement*. Springer, 2007, <https://doi.org/10.1007/978-1-4020-5302-3>.
30. Soulet, Guillaume, et al. “Black Sea ‘Lake’ Reservoir Age Evolution since the Last Glacial - Hydrologic and Climatic Implications.” *Earth and Planetary Science Letters*, vol. 308, no. 1-2, 2011, pp. 245-258, <https://doi.org/10.1016/j.epsl.2011.06.002>.
31. Micklin, Philip. “The Aral Sea Disaster.” *Annual Review of Earth and Planetary Sciences*, vol. 35, 2007, pp. 47-72, <https://doi.org/10.1146/annurev.earth.35.031306.140120>.
32. Micklin, Philip. “The Future Aral Sea: Hope and Despair.” *Environmental Earth Sciences*, vol. 75, no. 9, 2016, 844, <https://doi.org/10.1007/s12665-016-5614-5>.
33. Leblanc, Marc J., et al. “Evidence for Megalake Chad, North-Central Africa, during the Late Quaternary from Satellite Data.” *Palaeogeography, Palaeoclimatology, Palaeoecology*, vol. 230, no. 3-4, 2006, pp. 230-242, <https://doi.org/10.1016/j.palaeo.2005.07.016>.
34. Armitage, Simon J., et al. “West African Monsoon Dynamics Inferred from Abrupt Fluctuations of Lake Mega-Chad.” *Proceedings of the National Academy of Sciences*, vol. 112, no. 28, 2015, pp. 8543-8548, <https://doi.org/10.1073/pnas.1417655112>.
35. Coe, Michael T., and Jonathan A. Foley. “Human and Natural Impacts on the Water Resources of the Lake Chad Basin.” *Journal of Geophysical Research: Atmospheres*, vol. 106, no. D4, 2001, pp. 3349-3356, <https://doi.org/10.1029/2000JD900587>.
36. Pham-Duc, Binh, et al. “The Lake Chad Hydrology under Current Climate Change.” *Scientific Reports*, vol. 10, 2020, 5498, <https://doi.org/10.1038/s41598-020-62417-w>.
37. AghaKouchak, Amir, et al. “Aral Sea Syndrome Desiccates Lake Urmia: Call for Action.” *Journal of Great Lakes Research*, vol. 41, no. 1, 2015, pp. 307-311, <https://doi.org/10.1016/j.jglr.2014.12.007>.
38. Pekel, Jean-François, et al. “High-Resolution Mapping of Global Surface Water and Its Long-Term Changes.” *Nature*, vol. 540, no. 7633, 2016, pp. 418-422, <https://doi.org/10.1038/nature20584>.
39. Candiani, Vera S. *Dreaming of Dry Land: Environmental Transformation in Colonial Mexico City*. Stanford University Press, 2014, <https://doi.org/10.1515/9780804791076>.
40. Dong, Zhibao, et al. “Research Progress in China's Lop Nur.” *Earth-Science Reviews*, vol. 111, no. 1-2, 2012, pp. 142-153, <https://doi.org/10.1016/j.earscirev.2011.11.003>.
41. United States Geological Survey. *TopoView: Historical Topographic Maps*. USGS, <https://ngmdb.usgs.gov/topoview/>.
42. United States Geological Survey. *National Water Information System*. USGS, <https://waterdata.usgs.gov/>.
43. Lambeck, Kurt. “Shoreline Reconstructions for the Persian Gulf since the Last Glacial Maximum.” *Earth and Planetary Science Letters*, vol. 142, no. 1-2, 1996, pp. 43-57, <https://doi.org/10.1016/0012-821X(96)00069-6>.
44. Heyvaert, Vanessa M. A., and Cecile Baeteman. “Holocene Sedimentary Evolution and Palaeocoastlines of the Lower Khuzestan Plain (Southwest Iran).” *Marine Geology*, vol. 242, no. 1-3, 2007, pp. 83-108, <https://doi.org/10.1016/j.margeo.2007.01.008>.
45. Pournelle, Jennifer R. *Marshland of Cities: Deltaic Landscapes and the Evolution of Early Mesopotamian Civilization*. University of California, San Diego, 2003, <https://search.worldcat.org/search?q=Marshland+of+Cities+Pournelle>.
46. Stanley, Daniel Jean, and Andrew G. Warne. “Worldwide Initiation of Holocene Marine Deltas by Deceleration of Sea-Level Rise.” *Science*, vol. 265, no. 5169, 1994, pp. 228-231, <https://doi.org/10.1126/science.265.5169.228>.
47. Amorosi, Alessandro, et al. “Global Sea-Level Control on Local Parasequence Architecture from the Holocene Record of the Po Plain, Italy.” *Marine and Petroleum Geology*, vol. 87, 2017, pp. 99-111, <https://doi.org/10.1016/j.marpetgeo.2017.01.020>.
48. Herrmann, Albert. *Historical and Commercial Atlas of China*. Harvard University Press, 1935, <https://search.worldcat.org/search?q=Historical+and+Commercial+Atlas+of+China+Herrmann>.
49. Tan, Qixiang, editor. *The Historical Atlas of China*. 8 vols., SinoMaps Press, 1982-1987, <https://search.worldcat.org/search?q=Zhongguo+lishi+ditu+ji+Tan+Qixiang>.
50. Zong, Yongqiang, et al. “An Evolutionary Model for the Holocene Formation of the Pearl River Delta, China.” *The Holocene*, vol. 19, no. 1, 2009, pp. 129-142, <https://doi.org/10.1177/0959683608098957>.
51. Sinsakul, Sin. “Late Quaternary Geology of the Lower Central Plain, Thailand.” *Journal of Asian Earth Sciences*, vol. 18, no. 4, 2000, pp. 415-426, <https://doi.org/10.1016/S1367-9120(99)00075-9>.
52. Tanabe, Susumu, et al. “Stratigraphy and Holocene Evolution of the Mud-Dominated Chao Phraya Delta, Thailand.” *Quaternary Science Reviews*, vol. 22, no. 8-9, 2003, pp. 789-807, <https://doi.org/10.1016/S0277-3791(02)00242-1>.
53. Nguyen, Van Lap, et al. “Late Holocene Depositional Environments and Coastal Evolution of the Mekong River Delta, Southern Vietnam.” *Journal of Asian Earth Sciences*, vol. 18, no. 4, 2000, pp. 427-439, <https://doi.org/10.1016/S1367-9120(99)00076-0>.
54. Tamura, Toru, et al. “Initiation of the Mekong River Delta at 8 ka: Evidence from the Sedimentary Succession in the Cambodian Lowland.” *Quaternary Science Reviews*, vol. 28, no. 3-4, 2009, pp. 327-344, <https://doi.org/10.1016/j.quascirev.2008.10.010>.
55. Sengupta, Dhritiraj, et al. “Mapping 21st Century Global Coastal Land Reclamation.” *Earth's Future*, vol. 11, no. 2, 2023, e2022EF002927, <https://doi.org/10.1029/2022EF002927>.
56. Seasholes, Nancy S. *Gaining Ground: A History of Landmaking in Boston*. MIT Press, 2003, <https://mitpress.mit.edu/9780262194945/gaining-ground/>.
57. Zainal, Khadija, et al. “The Cumulative Impacts of Reclamation and Dredging on the Marine Ecology and Land-Use in the Kingdom of Bahrain.” *Marine Pollution Bulletin*, vol. 64, no. 7, 2012, pp. 1452-1458, <https://doi.org/10.1016/j.marpolbul.2012.04.004>.
58. Riding, Tim. “‘Making Bombay Island’: Land Reclamation and Geographical Conceptions of Bombay, 1661-1728.” *Journal of Historical Geography*, vol. 59, 2018, pp. 27-39, <https://doi.org/10.1016/j.jhg.2017.08.005>.
59. Rijkswaterstaat. *Zuiderzeewerken and Deltawerken*. Rijkswaterstaat, <https://www.rijkswaterstaat.nl>.
60. Port of Rotterdam Authority. *Maasvlakte 2*. Port of Rotterdam, <https://www.portofrotterdam.com>.
61. Geospatial Information Authority of Japan. *Maps and Geospatial Information*. GSI, <https://www.gsi.go.jp>.
62. Urban Redevelopment Authority. *Land Use and Reclamation Records*. URA Singapore, <https://www.ura.gov.sg>.
63. Lands Department. *Hong Kong Land Reclamation Records*. Government of the Hong Kong SAR, <https://www.landsd.gov.hk>.
64. Direcção dos Serviços de Cartografia e Cadastro. *Macau Land Area Records*. DSCC, <https://www.dscc.gov.mo>.
65. Korea Rural Community Corporation. *Reclamation Projects*. KRC, <https://www.ekr.or.kr>.
66. United States Geological Survey. *San Francisco Bay Historical Shoreline*. USGS, <https://www.usgs.gov>.
67. Suez Canal Authority. *Canal History*. SCA, <https://www.suezcanal.gov.eg>.
68. Panama Canal Authority. *Canal History*. ACP, <https://pancanal.com>.
69. Klein Goldewijk, Kees, et al. “Anthropogenic Land Use Estimates for the Holocene - HYDE 3.2.” *Earth System Science Data*, vol. 9, no. 2, 2017, pp. 927-953, <https://doi.org/10.5194/essd-9-927-2017>.
70. Schiavina, Marcello, et al. *GHS-POP R2023A - GHS Population Grid Multitemporal (1975-2030)*. European Commission, Joint Research Centre, 2023, <https://doi.org/10.2905/2FF68A52-5B5B-4A22-8F40-C41DA8332CFE>.
71. Kosarev, Aleksey N., et al. “Kara-Bogaz-Gol Bay: Physical and Chemical Evolution.” *Aquatic Geochemistry*, vol. 15, no. 1-2, 2009, pp. 223-236, <https://doi.org/10.1007/s10498-008-9054-z>.
72. Dumont, Henri J. “The Caspian Lake: History, Biota, Structure, and Function.” *Limnology and Oceanography*, vol. 43, no. 1, 1998, pp. 44-52, <https://doi.org/10.4319/lo.1998.43.1.0044>.
73. Kleiber, Hans Peter, and Frank Niessen. “Variations of Continental Discharge Pattern in Space and Time: Implications from the Laptev Sea Continental Margin, Arctic Siberia.” *International Journal of Earth Sciences*, vol. 89, no. 3, 2000, pp. 605-616, <https://doi.org/10.1007/s005310000130>.
74. Bauch, Henning A., et al. “Chronology of the Holocene Transgression at the North Siberian Margin.” *Global and Planetary Change*, vol. 31, no. 1-4, 2001, pp. 125-139, <https://doi.org/10.1016/S0921-8181(01)00116-3>.
75. Stein, Rüdiger. *Arctic Ocean Sediments: Processes, Proxies, and Paleoenvironment*. Developments in Marine Geology, vol. 2, Elsevier, 2008, <https://doi.org/10.1016/S1572-5480(08)X0001-6>.
76. Jakobsson, Martin, et al. “The International Bathymetric Chart of the Arctic Ocean Version 4.0.” *Scientific Data*, vol. 7, 2020, 176, <https://doi.org/10.1038/s41597-020-0520-9>.
77. Bahuguna, Anjali, et al. “IRS Views the Adams Bridge (Bridging India and Sri Lanka).” *Journal of the Indian Society of Remote Sensing*, vol. 31, no. 4, 2003, pp. 237-239, <https://doi.org/10.1007/BF03007343>.
78. Dandabathula, Giribabu, et al. “A High-Resolution Digital Bathymetric Elevation Model Derived from ICESat-2 for Adam's Bridge.” *Scientific Data*, vol. 11, no. 1, 2024, 705, <https://doi.org/10.1038/s41597-024-03550-3>.
79. Barnes, Richard, et al. “Priority-Flood: An Optimal Depression-Filling and Watershed-Labeling Algorithm for Digital Elevation Models.” *Computers & Geosciences*, vol. 62, 2014, pp. 117-127, <https://doi.org/10.1016/j.cageo.2013.04.024>.
80. Grant, K. M., et al. “Rapid Coupling between Ice Volume and Polar Temperature over the Past 150,000 Years.” *Nature*, vol. 491, no. 7426, 2012, pp. 744-747, <https://doi.org/10.1038/nature11593>.
81. Moody, John A., and Brent M. Troutman. “Characterization of the Spatial Variability of Channel Morphology.” *Earth Surface Processes and Landforms*, vol. 27, no. 12, 2002, pp. 1251-1266, <https://doi.org/10.1002/esp.403>.
82. Stein, Ruediger, et al. “Arctic (Palaeo) River Discharge and Environmental Change: Evidence from the Holocene Kara Sea Sedimentary Record.” *Quaternary Science Reviews*, vol. 23, no. 11-13, 2004, pp. 1485-1511, <https://doi.org/10.1016/j.quascirev.2003.12.004>.
83. Schwamborn, Georg, et al. “Late Quaternary Sedimentation History of the Lena Delta.” *Quaternary International*, vol. 89, no. 1, 2002, pp. 119-134, <https://doi.org/10.1016/S1040-6182(01)00084-2>.
84. Nobre, A. D., et al. “Height Above the Nearest Drainage – a Hydrologically Relevant New Terrain Model.” *Journal of Hydrology*, vol. 404, no. 1-2, 2011, pp. 13-29, <https://doi.org/10.1016/j.jhydrol.2011.03.051>.
  
</details>

<details>
  <summary>Trewartha (Expand to view sources).</summary>

1. Trewartha, Glenn T., and Lyle H. Horn. *An Introduction to Climate*. 5th ed., McGraw-Hill, 1980, <https://search.worldcat.org/search?q=An+Introduction+to+Climate+Trewartha+Horn+1980>.
2. Patton, Clyde P. “A Note on the Classification of Dry Climates in the Köppen System.” *The California Geographer*, vol. 3, 1962, pp. 105-112, <https://hdl.handle.net/10211.2/2248>.
3. Belda, Michal, et al. “Climate Classification Revisited: From Köppen to Trewartha.” *Climate Research*, vol. 59, no. 1, 2014, pp. 1-13, <https://doi.org/10.3354/cr01204>.
4. Liu, Zhengyu, et al. “Transient Simulation of Last Deglaciation with a New Mechanism for Bølling-Allerød Warming.” *Science*, vol. 325, no. 5938, 2009, pp. 310-314, <https://doi.org/10.1126/science.1171041>.
5. Otto-Bliesner, Bette, and Nan Rosenbloom. *Simulation of the Transient Climate of the Last 21,000 Years (TraCE-21ka)*. NSF National Center for Atmospheric Research, 2021, <https://doi.org/10.5065/CXB5-TV56>.
6. Fick, Stephen E., and Robert J. Hijmans. “WorldClim 2: New 1-km Spatial Resolution Climate Surfaces for Global Land Areas.” *International Journal of Climatology*, vol. 37, no. 12, 2017, pp. 4302-4315, <https://doi.org/10.1002/joc.5086>.
7. WorldClim. *WorldClim Version 2.1 Historical Climate Data (1970-2000, 5 Arc-Minute)*. University of California, Davis, 2020, <https://www.worldclim.org/data/worldclim21.html>.
8. Harris, Ian, et al. “Version 4 of the CRU TS Monthly High-Resolution Gridded Multivariate Climate Dataset.” *Scientific Data*, vol. 7, no. 1, 2020, 109, <https://doi.org/10.1038/s41597-020-0453-3>.
9. Climatic Research Unit, University of East Anglia. *CRU TS v4.10 Gridded Time-Series Dataset (1901-2025)*. University of East Anglia, 2026, <https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.10/>.
10. Tarasov, Lev, et al. “A Data-Calibrated Distribution of Deglacial Chronologies for the North American Ice Complex from Glaciological Modeling.” *Earth and Planetary Science Letters*, vol. 315-316, 2012, pp. 30-40, <https://doi.org/10.1016/j.epsl.2011.09.010>.
11. Briggs, Robert D., et al. “A Data-Constrained Large Ensemble Analysis of Antarctic Evolution since the Eemian.” *Quaternary Science Reviews*, vol. 103, 2014, pp. 91-115, <https://doi.org/10.1016/j.quascirev.2014.09.003>.
12. Ivanovic, Ruza F., et al. “Transient Climate Simulations of the Deglaciation 21-9 Thousand Years before Present (Version 1) - PMIP4 Core Experiment Design and Boundary Conditions.” *Geoscientific Model Development*, vol. 9, no. 7, 2016, pp. 2563-2587, <https://doi.org/10.5194/gmd-9-2563-2016>.
  
</details>
