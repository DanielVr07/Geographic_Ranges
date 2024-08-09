---
editor_options: 
  markdown: 
    wrap: sentence
---

# Estimating Geographic Ranges in Freshwater Species

R code for the analyses from "Comparing methods for estimating geographic ranges in freshwater fishes: several mirrors of the same reality" by Valencia-Rodríguez et al. 2024.
Freshwater Biology.
R code to replicate all data treatment and analyses from the paper

## File description:

**Data:** This folder contains the data from which the datasets used in all subsequent analyses were derived.

**R_Codes:** This folder contains three documents code.

**1. Extents of occurrence (EOO), is composed of five sections:**

1.1.
Optimal alpha value: This section calculates the optimal alpha value for the species automatically and adaptively, based on the density and distribution of the species' occurrence records.

1.2 Dynamic alpha: This section generates the species' geographical range using dynamic alpha values.
The polygons are overlapped on sub-basins and clipped to water bodies, resulting in restricted and unrestricted ranges to the bodies of water, as well as the sizes of the ranges.

1.3.
Static alpha: This section generates the species' geographical range using the same alpha value for all species (α = 6).
The polygons are overlapped on sub-basins and clipped to water bodies, resulting in restricted and unrestricted ranges to the bodies of water, as well as the sizes of the ranges.

1.4.
Convex hull: This section generates the species' geographical range using the convex hull method.
The polygons are also overlapped on sub-basins and clipped to water bodies, resulting in restricted and unrestricted ranges to the bodies of water, as well as the sizes of the ranges.

1.5 Expert maps: This section, expert maps available from the IUCN website are used to calculate range sizes, both unrestricted and restricted to water bodies.

**2. Species Distribution Models (SDMs) is composed of six sections:**

2.1.
Occurrence records: This section splits the occurrence data randomly for each species into training and validation sets.

2.2.
Variable selection: This section selects variables for each species using VIF and correlation analysis.

2.3.
Processing of the variables: In this section, the selected variables for each species are cropped based on their accessible area.

2.4.
Species Distribution Models: This section contains the parameters used to model the Curimatidae fish.

2.5.
Extraction of SDM Results: This section extracts the results of the SDM evaluations for each species and stores them in a data frame.

2.6.
Reclassification of the models: This section reclassifies the continuous models into binary models using the 10th percentile training presence value as the threshold.

**3. Analyses, is composed of three sections:**

3.1.
Comparison of range sizes between methods: This section compares the differences in estimated restricted geographic range sizes between methods.

3.2 Unrestricted vs. restricted range size within each polygon method: This section compares the range sizes obtained from constructing unrestricted polygons versus polygons restricted to freshwater bodies for each of the four polygon methods.

3.3.
Differences and similarities among species' range sizes: This section compares range sizes by evaluating pairs of methods at a time.

## References

García‐Andrade, A. B., Carvajal‐Quintero, J. D., Tedesco, P. A., & Villalobos, F.
(2021).
Evolutionary and environmental drivers of species richness in poeciliid fishes across the Americas.
Global Ecology and Biogeography, 30(6), 1245--1257.
<https://doi.org/10.1111/geb.13299>

Cobos, M. E., Townsend Peterson, A., Barve, N., & Osorio-Olvera, L.
(2019).
Kuenm: An R package for detailed development of ecological niche models using Maxent.
PeerJ, 2019(2), e6281.
<https://doi.org/10.7717/peerj.6281>

Barve, N., Barve, V., Jiménez-Valverde, A., Lira-Noriega, A., Maher, S. P., Peterson, A. T., Soberón, J., & Villalobos, F.
(2011).
The crucial role of the accessible area in ecological niche modeling and species distribution modeling.
Ecological Modelling, 222(11), 1810--1819.
<https://doi.org/10.1016/j.ecolmodel.2011.02.011>

Davis Rabosky, A. R., Cox, C. L., Rabosky, D. L., Title, P. O., Holmes, I. A., Feldman, A., & McGuire, J. A.
(2016).
Coral snakes predict the evolution of mimicry across New World snakes.
Nature Communications, 7(1), 11484.
<https://doi.org/10.1038/ncomms11484>

Domisch, S., Amatulli, G., & Jetz, W.
(2015).
Near-global freshwater-specific environmental variables for biodiversity analyses in 1 km resolution.
Scientific Data, 2(December), 1--13.
<https://doi.org/10.1038/sdata.2015.73>

Lehner, B., & Grill, G.
(2013).
Global river hydrography and network routing: Baseline data and new approaches to study the world's large river systems.
Hydrological Processes, 27(15), 2171--2186.
<https://doi.org/10.1002/hyp.9740>

## Contact

[davarod\@gmail.com](mailto:davarod@gmail.com){.email}, [fabricio.villalobos\@gmail.com](mailto:fabricio.villalobos@gmail.com){.email}, [octavio.rojas\@inecol.mx](mailto:octavio.rojas@inecol.mx){.email}
