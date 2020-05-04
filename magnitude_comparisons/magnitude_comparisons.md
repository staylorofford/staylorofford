# Magnitude Comparisons

## Project Summary

The magnitude comparison project aims to find the relationships between different types of earthquake magnitudes. As single earthquake catalogs are often limited by their purpose, this project includes functionality to match earthquakes between catalogs and compare their respective magnitudes.

## Motivations

Every magnitude type is limited in some way. The "simple" magnitudes like local or body wave magnitude - those derived directly from filtered seismogram amplitudes - are thought to saturate at high magnitudes and are often only approxmimate representations of the earthquake size, meaning they are only useful over a certain range of earthquake sizes and then have limitations in their accuracy. However, these simple magnitudes are calculated directly from the earthquake seismogram, meaning there is very little delay between seismogram observation and magnitude estimation with these magnitudes.  

The moment magnitude is by comparison a "slow" or "complex" magnitude: it requires modeling of the earthquake source using seismic waves that are often only recorded tens of minutes after the earthquake; it is this delay between seismogram observation and magnitude determination that makes this a relatively slow magnitude to compute. While faster variants of the moment magnitude exist, they are still produced minutes after the simple magnitudes would be for an earthquake. That said, the moment magnitude is the "truest" indication of an earthquake's size, meaning it is the desired magnitude in all cases. 

In earthquake response and tsunami threat assessment scenarios, the time taken to produce a moment magnitude is too long. Moment magnitude estimates are produced after the first, and sometimes the second, iteration of threat assessment. This can lead to delays in threat assessment if the scientist providing such assessment chooses to wait for the moment magnitude, or it can cause confusion and "derailing" of the response due to the moment magnitude being adopted part way through the response and forcing a "back-peddling" of the response if the threat assesment changes as a result. 

In reality, the second of the two cases presented above is more likely. When the simple magnitudes are available synchronously with the  earthquake location (i.e. within 5 minutes of the initial observation of the earthquake), the scientist performing threat assessment will be pressured - internally or externally - to use these magnitudes rather than wait for the moment magnitude. Sometimes the simple magnitude ends up being the same as the moment magnitude, a relieving situation. However, more often the simple magnitude and the moment magnitude differ, at times by up to half an order of magnitude. This mismatching, when occurring for earthquakes in the tsunamigenic magnitude range (say, M7+), can lead to gross under- or over-estimation of the tsunami potential of the earthquake.

These problems are non-trivial and can have serious consequences. The simple solution may appear to be to adjust the threat assessment process to wait for a moment magnitude, but if an earthquake of tsunamigenic size has occurred close to land (i.e. the tsunami will make landfall before the moment magnitude is calculated) then this is not an option. Another simple solution may appear to be to use the faster variant of the moment magnitude, but this suffers from the same problem, not to mention that many tsunami threat advisory agencies do not have the capability to produce moment magnitudes and are unable to gain such a capability. In reality, many such agencies rely on international institutions to provide moment magnitude estimates without any contract ensuring the availability of these services. Best efforts aside, this is not an ideal situation.

To summarise, there are "simple" magnitudes which are available at the same time as earthquake locations and there are "complex" magnitudes which, while much more accurate measures of earthquake size than simple magnitudes, are available at best many minutes after earthquake location. Often simple magnitudes are used following earthquakes of tsunamigenic size in tsunami threat assessment due to their availability and the immediate need of such threat assessment. When the moment magnitude, a "complex" magnitude, is available it is often used to update the threat assessment which, while improving assessment accuracy, can cause chaos in the response process. Adoption of a fast variant of the moment magnitude is possible, but may still be too slow when tsunamigenic earthquakes occur near land and may also require a risky dependence on international agencies to produce such magnitudes.

The current project aims to provide a bridge between the simple magnitudes and the moment magnitude so that estimates of the latter can be provided alongside the former. To do this, empirical relationships between the simple magnitudes and the moment magnitude are defined using data from 8 years of observation. Ultimately, the magnitude relationships are presented qualitatively in hopes of avoiding any over-dependence that could grow if functions describing the relationships are provided. It is the hope of the author that anyone using the relationships defined in this work does so knowing the limitations of the data used to define then. Moreover, the author hopes that such relationships will not be used in place of a moment magnitude; the moment magnitude should still be determined and used, these relationships are only intended to be used in the time between the earthquake location and the moment magnitude determination. However, it is of the greatest importance to the author that these relationships are used, and that their outputs help to provide accurate tsunami threat assessments and allow a smooth transition into use of the moment magnitude in the later stages of their associated responses. 

## Method

The magnitude comparisons project uses a collection of Python scripts to perform its calculations. The complete processing workflow, from catalog definition to relationship determination, follows.

### Catalog Definition

Using FDSN, a catalog is defined from an existing catalog using operator-defined paramters such as min/max latitude/longitude, depth, magnitude, and origin time. Up to two catalogs can be defined using different source catalogs with each being filtered by the same set of parameters.

### Comparison Magnitudes Definition

For each catalog the magnitude types to compare are defined. Each magnitude type will be compared against all others, including those in the same catalog.

### Matching Parameter Definition

Within the same catalog, earthquakes are matched by their unique ID.

Between different catalogs, matching follows the following process:

1. The euclidean distance between earthquake solutions in the distance-time plane is calculated by first calculating the euclidean distance between earthquake hypocentres in a cartesian coordinate system and then calculating the euclidean distance between this distance and the difference in the earthquake origin times. This gives a measure of how close in space and time two earthquake solutions (hypocentre and origin time) are.
1. For a given earthquake solution in the reference catalog, the euclidean distances to earthquake solutions in the comparison catalog are filtered by those whose relative origin time and relative spatial distance are within bounds set by the operator, e.g. only earthquake solutions whose relative origin time is within 100 seconds, and whose relative spatial distance is within 1000 km, of the reference earthquake solution are considered for matching. While not ultimately required, this filtering acts to reduce the computation time required for matching.
1. All earthquake solutions from the comparison catalog that meet the filters are then put through an inversion routine: using the earthquake solution from the reference catalog and the arrival times for the earthquake solution in the reference catalog, the RMS error between predicted arrival times from the comparison catalog earthquake solution and the arrival times from the reference catalog earthquake solution is calculated using the IASPEI91 travel times and the station locations of those stations with arrival times for the earthquake solution in the reference catalog.
1. The earthquake solution from the comparison catalog that gives the minimum RMS error is considered the match of the earthquake solution in the reference catalog, granted that the RMS error is below an operator-defined threshold, e.g. 5 seconds.

Once a matched earthquake solution is found, the magnitudes of that solution are added to those of the reference catalog earthquake solution. The result is a single catalog of earthquake solutions and all the magnitudes - when existing - of the desired types for each event in that catalog. If no comparison catalog is defined, the analysis dataset will consist only of the data from the reference catalog. Similarly, if no viable matches exist for earthquake solutions in the reference catalog then the analysis dataset will consist of only the data from the reference catalog.

### Analysis Dataset

The analysis dataset is derived from the matched earthquake solutions in up to two catalogs. It consists of earthquake unique IDs, hypocentre parameters, origin times (all from the reference catalog), and magnitude values for the desired magnitude types (from all catalogs), when they exist.

### Magnitude Comparison

Using the analysis dataset, the values of each combination of magnitude types is plotted in 2-dimensions both as a scatter plot of all value pairs and as an overlaid scatter plot of the average value of the reference catalog magnitude for each 0.1 unit bin of comparison catalog magnitude values (and the size of these points is proportional to the number of values in each bin). The line with intercept 0 and slope 1 (representing a perfect 1:1 relationship between the magnitudes) is then overlaid for reference. Lastly, a grid is overlaid to assist in determining the magnitude relationships.

The idea is to produce figures that can be used by a scientist to qualitatively determine the relationship between two magnitude types for a given value of one of those magnitude types. While it is possible to derive functions describing these relationships, this has not been included in the final iteration of this work for reasons explored in the project motivations.

Each plot is saved directly as a PNG file with axis labels noting the magnitude types being compared and a title noting the number of data points on the plot.

    ### Gutenberg-Richter Plots
    
    The author recognises that the datasets used in the magnitude comparions may be skewed by the distribution of values of each magnitude type. To mitigate the effect of such data featues Gutenberg-Richter style plots are included in the script output. These plots show the number of magnitude values in a given 0.1 unit bin. In most cases, the expectation is that 10x fewer earthquakes will exist for a 1 unit increase in magnitude value over the range in which the earthquake catalog is complete. The distribution of magnitude values often appears like the first half of a normal distribution at low values, then as a linearly decreasing line over the values for which the catalog is complete, and then as a line asymptotically decreasing towards 0 for high magnitudes. Deviations from this trend, such as bumps or bimodal distributions, suggest some kind of distortion in the data. While further exploration of such distortions is not considered as part of this work, knowledge of such distortions is important for when the scientist uses the magnitude comparison plots to determine magnitude relationships as they point to areas of the data that may be unrepresentative and as such should be could skew any interpretations. 

## Results

### Magnitude Comparisons Script

Presented below are the magnitude comparison plots produced by the `magnitude_comparisons.py` script. The script is available on the author's GitHub repository and contains in-line documentation for its use, features, and logic. The documentation is by no means complete, but should be adequate for an operator with some skill in programming. As much as the plots presented below are a result of this work, so is this script, and it is the hope of the author that this script can be used to reproduce the plots that follow in an improved format and with more data in the future.

### GeoNet Magnitudes versus Moment Magnitude

A unified moment magnitude dataset was created by merging moment magnitude values from the USGS earthquake catalog with those of the GeoNet moment tensor catalog. When two values existed for an event, the USGS moment magnitude was taken. In doing this, the author assumes that the two moment magnitude values are the same, an assumption that is supported by the work of John Ristau (REF) and which can be qualitatively verified from Figure 1 below. Each of the GeoNet magnitudes and their relationship with the moment magnitude are presented in the plots below with y-axis labels denoting which GeoNet magnitude is being compared. In addition to each magnitude comparions plot, a plot of the data distribution for each dataset is provided.



### GeoNet Magnitudes versus Moment Magnitude in SeisComP3 Regions

Often magnitude relationships vary at the regional scale with some areas producing moment magnitude estimates lower than those of the simple magnitudes and with other areas producing the opposite. To look at magnitude relationships in different regions the data is split by SeisComP3 region, this being given as a "description" in each event's entry in their respective earthquake catalog when queried via FDSN. In addition to each magnitude comparions plot, a plot of the data distribution for each dataset is provided.

## Discussion

- Apparent relationships
- Saturation
- Scale validity: low and high magnitudes
- Data distribution
- Take home message

## Conclusion

- Provide a reflection on the results and discussion in answer to the motivations

## References

- You know you need it: put it through the text and here
