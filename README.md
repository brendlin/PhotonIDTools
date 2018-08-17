PhotonIDTools - An Introduction to the Package
================

This package contains the following scripts:
 - **scripts/ValidateConfFiles.py**:
   - Use this script to validate whether Loose, Medium and Tight photon
conf files are subsets are one another, whenever validating any
new cut-based ID.
 - **scripts/PlotEfficiencies.py**:
   - This script can be used to plot cut-based efficiencies.
 - **scripts/PlotShowerShapeVariables.py**:
   - This script can be used to plot shower shape variables and the cut values applied to them. 
 - MC-to-MC AF2 scale factors (coming soon)


**PlotEfficiencies.py** - Description and Instructions
==================

### What is it
The **scripts/PlotEfficiencies.py** script is
used for evaluating (p<sub>T</sub>-dependent or independent) cut-based efficiencies,

### How to run it
The script is called in the following way:

```bash
PlotEfficiencies.py --menu1 Tight.conf --menu2 MyMenu.conf --names Tight,p_{T}-dependent --singlephotonsignal sp.root --radzsignal zy.root --outdir OutputDir
```

Where **sp.root** is a root file produced by the 
[HGamSinglePhotons](https://gitlab.cern.ch/ATLAS-EGamma/Software/PhotonID/HGamSinglePhotons "HGamSinglePhotons") package, and
Where **zy.root** is a root file produced by the 
[RadiativeZ](https://gitlab.cern.ch/ATLAS-EGamma/Software/PhotonID/RadiativeZ "RadiativeZ") package.

The background from jet-filtered samples (or others from the **HGamSinglePhotons** package)
can be plotted using the `--jetfiltered` option (we suggest running on the jet-filtered sample
separately, since the scale of the efficiency plot will be very different from the signal).

The script supports plotting up to four cut-based menus at a time (`--menu1` to `--menu3`),
with names given using the `--names` option (using comma-separated values,
e.g. `--names 'Tight,Loose,Menu 3 name,Menu 4 name'`).

### The output

The output of the script is a plot of the efficiencies as a function of
p<sub>T</sub> in each of the &eta; bins. Part of the goal of the script is
to reduce the number of output plots to a manageable number, to reduce "plot fatigue."
(The number of plots is reduced from 28 &rarr; 4.)

**PlotShowerShapeVariables.py** - Description and Instructions
==================

### What is it
The **scripts/PlotEfficiencies.py** script is
used to visualize shower shape variables and the rectangular cuts applied.

### How to run it
The script is called in the following way:

```
PlotShowerShapeVariables.py --tight Tight.conf --loose Medium.conf --names 'Tight ID,Medium ID' --singlephotonsignal SP.root  --outdir OutputDir
```

The script supports plotting up to four cut-based menus at a time
(The first two, `--tight` and `--loose`, are suggestively labeled,
and the second two `--menu3` and `--menu4`, with names given using the `--names` option
(using comma-separated values, e.g. `--names 'Tight,Loose,Menu 3 name,Menu 4 name'`).

Other options are available, including running simultaneously over other files:
 - **--radzsignal MyRadZSignalFile.root**
 - **--radzdata MyRadZDataFile.root**
 - **--singlephotonsignal MySPSignalFile.root**
 - **--singlephotondata MySPDataFile.root**
 - **--jetfiltered MyJFFile.root**
 - **--FixedCutLoose** (this is a boolean, to apply a FixedCutLoose isolation preselection)
 
### The output

The output of this script is a set of pdfs containing the plotted shower shapes
and the cuts applied to them. To reduce "plot fatigue", visually similar plots are grouped
to focus on the evolution of shower shapes vs p<sub>T</sub>. Further grouping reduces the
number of pdfs from 1260 to 36, by a factor of 35.
