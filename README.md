# Materials for People Need a Few Seconds to Be Random
This repository contains data and materials for People Need a Few Seconds to Be Random. 

The project is maintained by Lucas Castillo, please contact me at [lucas.castillo-marti@warwick.ac.uk](mailto:lucas.castillo-marti@warwick.ac.uk) or [castillo.lucas@protonmail.com](mailto:castillo.lucas@protonmail.com).

## Running the code
### Getting Started
The project ran on R version 4.3.3. It uses the `renv` for reproducible environments (see [here](https://rstudio.github.io/renv/articles/renv.html)). To get started: 
1. Open the `afstbr.Rproj` file with RStudio.
1. Rstudio will detect that it is a 'renv' project, and download and install `renv` for you. 
1. You can run `renv::restore()` in the console to install the packages needed. It will install the same versions we used. These will go in a separate folder (thus not affecting your other projects).

### Outline
Analyses we ran live in the main folder, and are numbered in the order they should be run (i.e., 1a, 1b, 2a, ...). They're clustered in three broad themes:
1. Data preparation and cleaning
1. Random forest training (simulations, etc.)
1. Applying forests to the data.

Some analysis files are separate to this main pipeline: 
- Analyses for Appendix A are named `AppA___.R` and can be run after `1a_exclusions.R`.
- Code to produce figures is named `P______.R` and will need for the main pipeline to have run. 
- The `O_demographics.R` file gets gender and age information for reporting (after `1a_exclusions.R`).

Other important folders: 
- `data/` (see below for license information)
- `manual_figures/` has the `.svg` we made for Figure 1 and the icons used to produce it (see below for license information)
- Some needed functions are in `src/` file. Note that `compute_all_rand_measures.R` is adapted from Angelike & Musch, 2024 (see below for license information).
- The `output/` folder stores:
    - `descriptive_analyses/`: results of the Appendix A analyses
    - `plots/` figures produced programatically (by `P___.R` files). 
- The `renv/` folder stores things needed by the `renv` package to function. 

## Cite this project  <!--bibtex -->
```
@{
    ...
}
```

## License 
This project is distributed under a [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/deed.en) license. However, note that some components present here (some data, some code) are not authored by us and thus might use different licenses.

Information about data sources and their licenses can be found [in the data folder](data/README.md). Information about icons used in Figure 1 can be found in the [manual_figures folder](manual_figures/icons/README.md). Finally, code in `src/compute_all_rand_measures.R` is adapted from [Angelike & Musch, 2024](https://doi.org/10.3758/s13428-024-02456-7) and distributed under a [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/deed.en) license.
