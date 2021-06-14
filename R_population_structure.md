## Population structure analysis R packages:

### LEA -- landscape and ecological association (Frichot and Francois, 2015, Methods in Ecology and Evolution; also see Francois, 2016 "Running Structure-like population genetic analyses with R" tutorial)
- performs pop structure analysis
    - different algorithm as STRUCTURE but claims similar estimates of ancestry (i.e., admixture) coefficients for out-crossing species and more accurate estimates in the presence of inbreeding
- number of clusters chosen by cross-entropy criterion (same as ADMIXTURE software; see tutorial)
- Same package can be used for genome scans of adative alleles and latent factor mixed modeling (LFMM) for ecological association analyses
- uses .lfm and .geno data formats and provides methods to convert from VCF

### pophelper (Francis, 2017, Mol Ecol Res; cited by Pradis et al. 2017 "integrated R ecosystem" paper)
- primarily used for visualization of results from other software (STRUCTURE, BAPS, TESS, etc)

### tess3r (Caye et al. 2017, Mol Ecol Res; extensions of LEA)
- incorporates geographic information to population structure analysis
- complementary to LEA package (with some overlapping functions)
- TESS appropriate when "the levels of ancestral population divergence are low"
- No HWE assumption
