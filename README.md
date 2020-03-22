<b> Citation</b>: Barnes, C. L., A. H. Beaudreau, M. W. Dorn, K. K. Holsman, and F. J. Mueter. 2020. Development of a predation index to assess trophic stability in the Gulf of Alaska. Ecological Applications. 

## Overview
This repository details the methods used to calculate indices of predation for Walleye Pollock (<i>Gadus chalcogrammus</i>) in the Gulf of Alaska (MT per year; 1990 to 2015). Pollock predators included: Arrowtooth Flounder (<i>Atheresthes stomias</i>), Pacific Cod (<i>Gadus macrocephalus</i>), Pacific Halibut (<i>Hippoglossus stenolepis</i>), Sablefish (<i>Anoplopoma fimbria</i>), and Walleye Pollock conspecifics. We used predation indices to estimate synchrony in pollock consumption and make inferences about trophic stability among demersal fishes in the Gulf of Alaska.

## File Structure
Input data (survey and food habits), shapefiles, and model outputs can all be found in Folder 1 ('1_Data' folder). Folders 2 through 5 contain species-specific analyses for each component of the predation index (resulting estimates are stored in the data folder). Folder 6 contains predation indices (all predators combined) for each of the spatial scales of interest: basin, the area encompassed by the stock assessment for Gulf of Alaska pollock, subregion, and statistical area. Script files in Folder 6 also include variance ratio calculations, which enabled estimates of synchrony and portfolio effects. Specific analyses that resulted in publication tables and figures are noted throughout. 

## Data Sources
<b>Total Predator Biomass</b>: Total biomass estimates were obtained from the most recent stock assessment for each groundfish predator (Barbeaux et al. 2017, Dorn et al. 2017, Hanselman et al. 2017, Spies et al. 2017, Stewart and Hicks 2017). Coast-wide estimates for Pacific Halibut were adjusted to reflect biomass in the Gulf of Alaska.

<b>Relative Predator Densities</b>: Bottom trawl survey data (all groundfish predators; 1990 to 2017) were collected by the Resource Assessment and Conservation Engineering (RACE) Division of the Alaska Fisheries Science Center (AFSC, NOAA) and are publicly accessible at https://www.afsc.noaa.gov/RACE/groundfish/survey_data/data.htm. See von Szalay et al. (2016) for information about bottom trawl survey design and data collection methods. Setline survey data (Pacific Halibut; 1998 to 2017) were collected by the International Pacific Halibut Commission and are publicly available at: https://iphc.int/data/fiss-data-query. For setline survey methods, see Clark and Hare (2006). Longline survey data (Sablefish; 1990 to 2017) were collected by the AFSC's Auke Bay Laboratories and can be found at https://www.afsc.noaa.gov/maps/longline/Map.php. See Sigler and Zenger (1989) for methods descriptions of the longline survey. 

<b>Mean Annual Rations and Age-specific Proportions of Pollock Consumed</b>: Food habits data (all groundfish predators; 1990 to 2015) were provided by the AFSC's Resource Ecology and Ecosystem Modeling (REEM) Program and are publicly accessible at: https://access.afsc.noaa.gov/REEM/WebDietData/DietDataIntro.php. For food habits data collection and processing methods, see Livingston et al. (2017).  

## Financial and Logistical Support
This project was funded by the Pollock Conservation Cooperative Research Center (G00009488) and the Rasmuson Fisheries Research Center associated with the University of Alaska Fairbanks. An anonymous donor supplied additional funds via the Northern Gulf of Alaska Applied Research Award. The University of Alaska (Juneau Fisheries Division and Southeast Sitka Campus) provided facilities and additional support. 

## Acknowledgments
We appreciate assistance with data acquisition and processing from Kerim Aydin, Steve Barbeaux, Troy Buckley, Dana Hanselman, Tom Kong, Geoff Lang, Wayne Palsson, and Ian Stewart. Jordan Watson and Lorenzo Ciannelli provided guidance on the initial development of spatial models. Mary Hunsicker and two anonymous reviewers provided valuable comments to improve upon the analyses detailed here. <br>

The authors would like to acknowledge Terry Quinn for offering his insight and expertise during earlier stages of this project. We have dedicated this work to him.

## References 

### Stock Assessments
&#8212; Barbeaux, S., K. Aydin, B. Fissel, K. Holsman, and W. Palsson. 2017. Assessment of the Pacific cod stock in the Gulf of Alaska. North Pacific Fishery Management Council Gulf of Alaska SAFE Report 189–332. <br>
&#8212; Dorn, M., K. Aydin, B. Fissel, D. Jones, A. McCarthy, W. Palsson, and K. Spalinger K. 2017. Assessment of the Walleye Pollock stock in the Gulf of Alaska. North Pacific Fishery Management Council Gulf of Alaska SAFE Report 47–182. <br>
&#8212; Hanselman, D. H., C. J. Rodgveller, C. R. Lunsford, and K. H. Fenske. 2017. Assessment of the Sablefish stock in Alaska. North Pacific Fishery Management Council Bering Sea, Aleutian Islands, and Gulf of Alaska SAFE Report 327–502. <br>
&#8212; Spies, I., K. Aydin, J. N. Ianelli, and W. Palsson. 2017. Assessment of the Arrowtooth Flounder stock in the Gulf of Alaska. North Pacific Fishery Management Council Gulf of Alaska SAFE Report 749–846. <br>
&#8212; Stewart, I., and A. Hicks. 2017. Assessment of the Pacific halibut (Hippoglossus stenolepis) stock at the end of 2017. International Pacific Halibut Commission IPHC-2018-AM094-10. <br>
### Survey and Food Habits Data
&#8212; Clark, W. G., and S. R. Hare. 2006. Assessment and management of Pacific halibut: data, methods, and policy. IPHC Scientific Report 83. <br> 
&#8212; Livingston, P. A., K. Aydin, T. W. Buckley, G. M. Lang, M-S. Yang, and B. S. Miller. 2017. Quantifying food web interactions in the North Pacific – a data-based approach. Environmental Biology of Fishes. 100(4):443–470. <br>
&#8212; Sigler, M. F., and H. H. Zenger, Jr. 1989. Assessment of Gulf of Alaska Sablefish and other groundfish based on the domestic longline survey, 1987. NOAA Technical Memorandum NMFS-AFSC Report 169. <br>
&#8212; von Szalay, P. G., and N. W. Raring. 2016. Data report: 2015 Gulf of Alaska bottom trawl survey. Seattle, WA. NOAA Technical Memorandum NMFS-AFSC-325. <br>
### Species Distribution Modeling
&#8212; Barnes, C. L., A. H. Beaudreau, M. E. Hunsicker, and L. Ciannelli (2018). Assessing the potential for competition between Pacific Halibut (Hippoglossus stenolepis) and Arrowtooth Flounder (Atheresthes stomias) in the Gulf of Alaska. PLoS ONE 13(12):e0209402. <br>
&#8212; Hunsicker, M. E., L. Ciannelli, K. M. Bailey, S. Zador, and L. Stige. 2013. Climate and demography dictate the strength of predator-prey overlap in a subarctic marine ecosystem. PLoS ONE 8(6):e66025. <br>
&#8212; Shelton, A. O., M. E. Hunsicker, E. J. Ward, B. E. Feist, R. Blake, C. L. Ward, et al. 2017. Spatio-temporal models reveal subtle changes to demersal communities following the Exxon Valdez oil spill. ICES Journal of Marince Science doi:10.1093/icesjms/fsx079. <br>
### Bioenergetics
&#8212; Armstrong JB and Schindler DE. 2011. Excess digestive capacity in predators reflects a life of feast and famine. Nature. 476:84–87. <br>
&#8212; Beaudreau, A. H., and T. E. Essington. 2009. Development of a new field-based approach for estimating consumption rates of fishes and comparison with a bioenergetics model for lingcod (Ophiodon elongatus). Canadian Journal of Fisheries and Aquatic Sciences 66:565−578. <br>
&#8212; Harvey, C. J. 2009. Effects of temperature change on demersal fisheries in the California Current: a bioenergetics approach. Canadian Journal of Fisheries and Aquatic Sciences 66:1449–1461. <br>
&#8212; Holsman, K. K., and K. Aydin. 2015. Comparative methods for evaluating climate change impacts on the foraging ecology of Alaskan groundfish. Marine Ecology Progress Series 521:217–235. <br>
&#8212; Holsman, K. K., K. Aydin, J. Sullivan, T. Hurst, and G. Kruse. 2019. Climate effects and bottom-up controls on growth and size-at-age of Pacific halibut (Hippoglossus stenolepis) in Alaska (USA). Fisheries Oceanography 28:345–358. <br>
### Miscellaneous
&#8212; Brodziak, J. 2012. Fitting length-weight relationships with linear regression using the log-transformed allometric model with bias-correction. NOAA Technical Memorandum PIFSC-H-12-03. <br>
&#8212; Chipps, S. R., and J. E. Garvey. 2007. Assessment of diets and feeding patterns. In: Analysis and interpretation of freshwater fisheries data. C. S. Guy and M. L. Brown, editors. Bethesda, MD. American Fisheries Society 473–514. 
