# Househunter

## Project Summary

The househunter project aims to find those properties in Wellington city that meet an operator-defined set of criteria. Due to the varying terrain of Wellington city and the various natural hazards that threaten any property in the area, such an aim is non-trivial. The househunter analytics code - the code which returns a shortlist of suitable properties from the operator-defined suitability criteria - combines many spatial datasets describing property features such as elevation, distance to coast, and sunlight. For each dataset, the operator defines the minimum value describing suitability for that dataset. Only those datasets which the operator defines a minimum value for are used in the property search. Output is in the form of a list of properties meeting the suitability criteria. Optionally a map can also be produced showing the location of the suitable properties.

## Motivations

Wellington city is a desirable place to live, but many of the properties are either dark or at risk of natural disasters. When looking to buy property for owner occupation, one desires to know more about the property than what can be gleamed from the documentation and a passing visit. The locations of suitable properties - those meeting a set of criteria such as average sunlight hours, elevation, etc. - often do not cluster, making the problem of finding a suitable property in the city non-trivial. The topography of the city cause the property features defining suitability, such as sunlight, to vary greatly between adjacent properties. Due to this, generalisations about any area in the city hold loosely at best, and to have any hope of accuracy in suitability estimates the quality of any property must be assessed independently of those nearby. The househunter project began so that this problem could be addressed. The goal of the project is to answer questions like, "Is this house going to get good sun? If it does, will it be safe from a tsunami?" Often times one cannot take for granted that such questions will be true, and when one is looking for a place to live for their lifetime such knowledge is of a high value. The motivation of the project is thus: that one may know, as derived from objective data, which properties in Wellington city meet one's definition of suitability. From this output one can then make decisions based on time-varying things, like personal finances, available properties, and so on.

## Suitability Criteria

The suitability criteria considered for the househunter project were defined exclusively from those of interest ot the project leader. They are presented as follows, with the property feature each describes included in parenthesis beside it:

- Average insolation (sunlight)
- Wind zone (shelter)
- Flooding risk (flooding risk)
- Ground shaking hazard (earthquake risk)
- Proximity to fault (earthquake risk)
- Slope (landslide risk)
- Elevation (sea level rise risk)
- Distance to coast (tsunami risk)
- Liquifaction potential (earthquake risk)

Each operator must define which suitability criteria they want to filter properties by, and above which value they consider a suitability criteria met. 

### Datasets

There is a dataset for each suitability criteria. Each dataset is derived from the combination of a non-property specific dataset of the Wellington city area describing each suitability criteria and a database of property tiles. To assign a single value to each property tile for each suitability criteria the minimum value of the suitability criteria within each property tile is used. Some suitability criteria datasets are imported without adjustment from the data source, while others are derived from such datasets.

#### Input Datasets

The input datasets used, presented with their source in parenthesis, in the order corresponding to the respective suitability criteria, are:

- *Average insolation (derived from LINZ 1m DSM)
- Wind Zones (WCC)
- Potential Flooding Hazard Area (WCC)
- Hazard Ground Shaking Area (WCC)
- Hazard Fault Line Area (WCC)
- *Slope (derived from LINZ 1m DEM)
- 1m DEM (LINZ)
- *Distance to coast (derived from LINZ NZ contours topo 150k)
- Liquefaction Potential (WCC)
- NZ property tiles (LINZ)

\* These datasets were derived. Details of how each dataset was derived is given in Appendix 1.

As all of these datasets are provided as, or derived from, government supplied datasets that come at no cost the project lead would like to express their gratefulness to the collectors, maintainers and providers of these datasets.

## Appendices

### Appendix 1: Dataset Derivation

Three datasets used in the househunter analysis code were derived as a part of the househunter project. Details on the derivation of each of these datasets is presented in this appendix.

#### Average Insolation

The input dataset used is the 1 m Digital Surface Model provided by LINZ and cut to the Wellington city area. The Digital Surface Model (DSM) was chosen over the Digital Elevation Model (DEM) due to its more accurate representation. While the DEM better captures topography, many of features which block sunlight in Wellington city are not topographic, e.g. houses, trees. Such features are not captured by the DEM, but are fairly reflected in the DSM. That said, such features in the DSM appear as the surface in sunlight models, which causes calculated insolation on some property tiles to be erronerously high. Ideally a combination of the DEM (to give the property tile surface) and the DSM (to give all surrounding obstructions to sunlight) would be used. In lieu of this, the DSM was used with the caveat that, when deriving average sunlight values for a property tile from DSM-derived sunlight models, a spatial averaging over the property tile may be required to mitigate the effect of surface location errors. This issue will be referred to hereafter as the "surface location problem".

To calculate sunlight models, GRASS GIS was used. First, aspect, slope, and horizon grids were calculated from the DSM. Horizon grids were calculated at 5 degree (20 minute) intervals about the horizontal plane. The calculation of these grids took 1.5 weeks using one processor of the project computer.

Once the aspect, slope, and horizon grids were complete, the insolation grids were calculated. These were made using the default parameters of `r.sun` and the pre-calculated grids in GRASS GIS. Due to the high processing load the code was run in low memory mode and the processing was split over 10 segments of the DSM for each day of processing. Using all four processors of the project computer the calculation took 1.5 weeks. Only the direct solar insolation was produced for each day of year processed. As the insolation on any DSM cell is sinusoidal in time and reflects about the shortest/longest day, only half the days of the year were processed to save disk space. This solar modeling was initially attempted without pre-computation of the aspect, slope, and horizon grids, but as this was estimated to take 3.5 years of computation time the pre-computation approach was favoured.

One insolation grid was produced for each day of the year. Grid values are the sum of insolation values over the day, these being calculated once every 20 minutes over the course of each day.

Due to the surface location problem, a simple spatial average over the insolation values on each property tile could not be used to calculate average sunlight values. Instead, the spatial average for each property tile would be calculated only over the subset of grid cells that could be considered representative of the property tile surface. It is assumed that those grid cells are the ones for which the elevation difference between the DSM and DEM on the tile is less than 1 metre. This threshold was set arbitrarily from inspection of the DSM-DEM difference grid and is assumed to represent the error introduced in the calculation of the DEM (which is derived from the DSM) and the sampling error introduced in the creation of the DSM and propagated forward into the DEM (the aggregation of elevation points from LIDAR into 1 m grid cells).

The average insolation for each grid cell was calculated by summing the net insolation value in that cell over each day of the year and then dividing by the number of days summed over. The average insolation value for a property tile was then calculated by summing the average insolation values in each grid cell within that property tile and dividing by the number of grid cells summed over. Only those grid cells considered representative of the property tile surface were summed over in the averaging.

#### Slope

#### Distance to Coast