PhotonIDTools
================

This package contains the following scripts:
 - `scripts/ValidateConfFiles.py`:
   - Use this script to validate whether Loose, Medium and Tight photon
conf files are subsets are one another, whenever validating any
new cut-based ID.
 - MC-to-MC AF2 scale factors (coming soon)
 
How to check out this package:
================
 ```
 git clone ssh://git@gitlab.cern.ch:7999/ATLAS-EGamma/Software/PhotonID/PhotonIDTools.git
 ```
 
How to add a file:
================
  - Put the file in the intended directory (e.g. `scripts`)
  - Do:
  ```
  git add scripts/myScript.C
  git commit -m "Adding myScript.C script"
  git push
  ```
  
**ScaleFactorPlots.py** - plot cosmetics for the Rel 21.2 Egamma Paper
================

I am adding just a brief description here of this script, which formats the scale factor plots for
the Rel 21.2 Egamma paper. The files themselves are still hard-coded -- if it becomes useful to
make this script more general, then it can be done in the future.

The usage is described below:
  ```
  ScaleFactorPlots.py --eta 0
  ScaleFactorPlots.py --eta 3 --converted
  ```
