Scripts to support the analysis of tunover for the USGS Climate Change and Biodiversity report.

No original data was generated for this project. All data is publicly available from:

BioTime: https://biotime.st-andrews.ac.uk/download.php ; https://doi.org/10.1111/geb.70003 (Note version 1 was used)

Marine Ecoregions of the World: https://databasin.org/datasets/3b6b12e7bcca419990c9081c0af254a2/ ; Spalding MD, Fox HE, Allen GR, Davidson N, Ferdaña ZA, Finlayson M, Halpern BS, Jorge MA, Lombana A, Lourie SA, Martin KD, McManus E, Molnar J, Recchia CA, Robertson J (2007) Marine Ecoregions of the World: a bioregionalization of coast and shelf areas. BioScience 57: 573-583.

Terrestrial Ecoregions: https://databasin.org/datasets/68635d7c77f1475f9b6c1d1dbe0a4c4c/ ; Olson, D.M., E. Dinerstein, E.D. Wikramanayake, N.D. Burgess, G.V.N. Powell, E.C. Underwood, J.A. D'Amico, I. Itoua, H.E. Strand, J.C. Morrison, C.J. Loucks, T.F. Allnutt, T.H. Ricketts, Y. Kura, J.F. Lamoreux, W.W. Wettengel, P. Hedao, and K.R. Kassem. Terrestrial Ecoregions of the World: A New Map of Life on Earth (PDF, 1.1M) BioScience 51:933-938.

Analyses on Google Earth Engine will require upload of shapefiles to generate distances to ecoregion boundaries.

Scripts should be run in the following order with the platform in ():
  1- Biodiversity Assessment Filtering BioTIME Data and Grids (R)
  2- Biodiversity Assessment Modeling Temperature and Precipitation Trends (R)
  3- Biodiversity Assessment Extract Cell Data (Google Earth Engine Code Editor)
  4- Run _targets file using targets::tar_make() (R)
  5- Analysis of turnover (R)

