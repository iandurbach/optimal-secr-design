# optimal-secr-design

Code for choosing detector locations for SCR surveys that maximize the approximate precision of density estimators. Detectors are placed to maximize whichever is the smaller of E(n), the number of animals detected (first captures), and E(r), the number of recaptures (total detections less first captures), based on the approximation in Efford and Boulanger (2019): CV(Dhat) ~ 1/sqrt[min{E(n),E(r)}]. The resulting designs are called min(n,r) designs. The code and output in this repo accompany the paper 

- Durbach, I., Borchers, D., Sutherland, C., & Sharma, K. Fast, flexible alternatives to regular grid designs for spatial capture-recapture. To appear in Methods in Ecology and Evolution.

Functions for constructing min(n,r) survey designs have been developed by combining two core design functions *scrdesignGA* and *scrdesignOF* from the R package *oSCR* with the *Enrm* function from the *secrdesign* package, which provides fast C implementation of the E(n) and E(r) calculations. Calculations for uniform and spatially-varying density can be handled by *secrdesign*`s *Enrm* function directly. For spatially-varying detection, we provide a modified version of this function -- this requires compiling the underlying C code. The *scrdesignGA* function in *oSCR* is essentially a wrapper function around the *kofnGA* function in package *kofnGA*, which implements a genetic algorithm. 

Our design functions, renamed to *scrdesignGAenr* and *scrdesignOFenrm*, are available in the **oSCR** folder in this repository. The functions require the *secrdesign* package, primarily because the calculations of E(n) and E(r) use *secr* mask objects. Function *scrdesignOFenrm* is called by *scrdesignGAenrm* and requires no input from the user -- it has simply been updated to calculate the number of first captures and recaptures respectively. Function *scrdesignGAenrm* requires a number of inputs, these are described in the supplementary material of the paper above, or in the comments to the *scrdesignGAenrm.R* file in the *oSCR* folder.

The only files in this repo that are strictly required to generate min(n,r) designs are *scrdesignGAenrm.R* and *scrdesignOFenrm.R* in the **oSCR** directory. To use these, download the files into your project directory and source them with:

```
source(\my_project_dir\scrdesignGAenrm.R)
source(\my_project_dir\scrdesignOFenrm.R)
```

Required packages are *dplyr*, *sf*, *secrdesign*, and *oSCR*. For a simple example generating min(n,r) designs for uniform and non-uniform activity centre densities see *example-uniformD.R* and *example-spatialD.R* respectively. 

Running designs with spatially-varying *lambda_0* requires compilation of additional C++ code performing the E(n) and E(r) calculations. To do this, browse to your **oSCR** folder and, at the command line/terminal, enter `R CMD SHLIB mysecrdesign.c`. To source the compiled code within R, use `dyn.load("oSCR/mysecrdesign.so")`, as shown in the examples available in the repository. An example is given in *example-spatial-lambda0.R*. 

All other files reproduce the designs in our paper and results derived from these. Intermediate output objects used in the plots and tables in the paper are saved in the *output* folder.

- Efford, M., & Boulanger, J. (2019). Fast evaluation of study designs for spatially explicit capture-recapture. Methods in Ecology and Evolution, 10, 1529–1535.
- Sutherland, C., Royle, A., & Linden, D. (2018). oscr: Multi-session sex-structured spatial capture-recapture models [Computer software manual]. (R package version 0.42.0)
- Efford, M. G. (2019). secrdesign: Sampling Design for Spatially Explicit Capture-Recapture. R package version 2.5.7.
- Wolters, M.A (2015). A Genetic Algorithm for Selection of Fixed-Size Subsets with Application to Design Problems. Journal of Statistical Software, 68(1), 1-18. 