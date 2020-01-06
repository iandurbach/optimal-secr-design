# optimal-secr-design
Approximately optimal survey design for spatial capture-recapture

Code for choosing detector locations for SCR surveys that maximize the precision of density estimators. Detectors are placed to maximize whichever is the smaller of E(n), the number of animals detected (first captures), and E(r), the number of recaptures (total detections less first captures), based on the approximation in Efford and Boulanger (2019): CV(Dhat) ~ 1/sqrt[min{E(n),E(r)}].

Two functions *Qfn* (in *\oSCR\Qfn_mcvd.R*) and *SCRdesign* (in *\oSCR\SCRdesign_mcvd.R*) overwrite functions of the same name in the R package oSCR (S), modifying the objective function to maximize precision and providing some additional functionality for non-uniform activity centre densities and non-uniform space use. These are the only files that are strictly required to generate approximately optimal designs. To use these, download the files into your project directory and overwrite the equivalent oSCR functions with 

```
library(oSCR)
source(\my_project_dir\Qfn_mcvd.R)
source(\my_project_dir\SCRdesign_mcvd.R)
```

For a simple example generating approximately optimal designs for uniform and non-uniform activity centre densities see *example-uniformD.R* and *example-nonuniformD.R* respectively. 

All other files reproduce the designs in the paper "Approximately optimal survey design for spatial capture-recapture" (under review) and results derived from these. Intermediate output objects used in the plots and tables in the paper are saved in the *output* folder.

- Efford, M., & Boulanger, J. (2019). Fast evaluation of study designs for spatially explicit capture-recapture. Methods in Ecology and Evolution, 10, 1529–1535.
- Sutherland, C., Royle, A., & Linden, D. (2018). oscr: Multi-session sex-structured spatial capture-recapture models [Computer software manual]. (R package version 0.42.0)
