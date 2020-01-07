# Magnitude Comparisons

## Project Summary

The magnitude comparison project aims to find the relationships between different types of earthquake magnitudes. As single earthquake catalogues are often limited by their purpose, this project includes functionality to match earthquakes between catalogues and compare their respective magnitudes.

## Motivations

Every magnitude type is limited in some way. The "simple" magnitudes like local or body wave magnitude - those derived directly from filtered waveform amplitudes - saturate at high magnitudes and are often only roughly representative of the earthquake moment magnitude, this being the "truest" indication of an earthquake's size. However, the moment magnitude, while accurate, takes up to 20 minutes to calculate after an earthquake, delaying earthquake response and tsunami threat assessment. In comparison, "simple" magnitudes can be calculated as soon as the waveform is recorded at a station, giving timely indication of an earthquake's size. The magnitude comparisons project is thus motivated to find the bridge between the simple, fast magnitudes and the accurate moment magnitude so that the former can be readily translated into the latter during event response situations while the moment magnitude is being calculated. As any earthquake catalogue is limited, either by its purpose or simply by the region and time period it captures, the magnitude comparisons project requires that two earthquake catalogues can be "combined" so that the magnitudes in one can be compared to those in the other to create a more complete dataset. The example of this which drove the project development is that of combining a local earthquake catalogue containing only simple magnitudes with the global moment tensor catalogue so that the relationship between the magnitudes in the local catalogue and the moment magnitude could be determined. The benefit of having a large combined dataset is that magnitude saturation - as well as the magnitude relationships at non-saturated values - can be studied, allowing the suitability of any simple earthquake magnitude to be known for an earthquake of any given size.

## Method

The magnitude comparisons project uses a collection of Python code to perform its calculations. The complete processing workflow, from catalogue definition to relationship determination, follows.

### Catalogue Definition

Using FDSN, a catalogue is defined from an existing catalogue using operator-defined paramters such as min/max latitude/longitude, depth, magnitude, and origin time. Up to two catalogues can be defined using different source catalogues with each being filtered by the same set of parameters.

### Comparison Magnitudes Definition

For each catalogue the magnitude types to compare are defined. Each magnitude type specified will be compared against all others, including those within the same catalogue.

### Matching Parameter Definition

Currently (prior to the introduction of an inversion routine to test earthquake solutions from different catalogues for consistency), earthquake solution matching between catalogues is performed by minimising the euclidean distance between earthquake solutions in the relative distance-time plane. Filters on distance and time (these values are relative to the earthquake solution in the reference catalogue, so it is the difference in origin time and hypocentre between the reference catalogue earthquake solution and the solutions of the comparison catalogue earthquakes that is being minimised) are defined to exclude earthquakes solutions that are - as the operator sees it - too far in time and distance to be possible representations of the same earthquake. Optionally the operator can have the data in the relative distance-time plane plotted for each earthquake solution in the reference catalogue to inspect the data distribution and aid in matching parameter definition. Within these parameters, the comparison catalogue earthquake solution with the lowest euclidean distance is considered to match the reference catalogue earthquake solution, i.e. they are both representations of the same earthquake.

In future the euclidean distance minimisation will be replaced by an inversion routine to determine earthquake solution matches across catalogues. This will work by iterating through comparison catalogue earthquake solutions (sorted by the relative distance-time euclidean distance) for each reference catalogue earthquake solution and calculating the RMS travel time residual for each comparison catalogue earthquake solution using an earthquake location routine that is either the same as that used to define the reference catalogue earthquake solutions or one that is used to relocate the reference catalogue earthquakes prior to the solution matching routine. The comparison catalogue earthquake solution with the lowest RMS residual for each reference catalogue earthquake solution is chosen as the matched solution. As with the existing method, an operator-defined parameter of maximum RMS residual will require definition to ensure only earthquake solutions for the same earthquake are ever matched.

Once a matched earthquake solution is found, the magnitudes of that solution are added to those of the reference catalogue earthquake solution. The result is a single extended catalogue of earthquake solutions. If no comparison catalogue is defined, the analysis dataset will consist only of the data from the reference catalogue. Similarly, if no viable matches exist for earthquake solutions in the reference catalogue then the analysis dataset will consist of only the data from the reference catalogue.

### Analysis Dataset

As described in the definitions above, the analysis dataset is derived from the matched earthquake solutions in up to two catalogues. It consists of earthquake identifiers, hypocentres, origin times, and magnitude values for the desired magnitude types, when they can be provided.

### Magnitude Comparison

Using the analysis dataset, each combination of magnitude types is respectively plotted in 2-dimensions and an orthogonal projection is calculated between the data in each case. The line defined by the orthogonal projection is overlain on the plot and a c=0, m=1 line is included also for reference. Each plot is saved directly as a PNG file.

### Earthquake Probabilities

In addition to the magnitude comparison functionality, the current version of the magnitude comparison code can calculate earthquake probabilities using data from the analysis dataset and. The likelihood of an earthquake of the type in the analysis dataset is calculated by assuming that the earthquakes follow a Poisson distribution, i.e. that earthquakes occur randomly and the occurrence of one does not influence the occurrence of any other. Additional parameters are required to use this functionality: the time period to calculate earthquake likelihood within, e.g. 1 week, and the minimum number of earthquakes occurring within that time period to calculate the likelihood for, e.g. 1 (so here the question is: what is the likelihood of at least 1 earthquake of the type captured in the analysis dataset occurring within 1 week?). The earthquake likelihood is given for each magnitude type in the analysis dataset and for the largest possible magnitude of any of the types in the analysis dataset.

Optionally the operator can have the earthquake probabilites calculated and plotted over time by specifying a desired number of time bins and overlap between the bins. The plot produced will indicate whether earthquake likelihood is time-varying and how it changes over time. The earthquake likelihood is here calculated using the largest possible magnitude for each earthquake solution in the analysis data. 

These outputs are produced using a highly simplified and generalised method of earthquake probability calculation and are not intended to be used in any capacity beyond informing the operator and acting as a tool for quality-checking the magnitude comparisons (for example, the operator may ask: is it fair to include certain time periods if the earthquake likelihood varies over that time period? Specifically, are the magnitude relationships influenced by phenomena such as mainshock-aftershock sequences or swarms? Or by changes in the instrument network or detection system used to define one or both of the catalogues?). 