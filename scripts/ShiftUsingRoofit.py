#!/usr/bin/env python

import os,sys
import ROOT
from array import array

# Add base path
the_path = ('/').join(os.path.abspath(__file__).split('/')[:-2]) 

python_path = '%s/python'%(the_path)
print 'Adding %s to PYTHONPATH.'%(python_path)
sys.path.append(python_path)

genericUtils_path = '%s/genericUtils'%(the_path)
print 'Adding %s to PYTHONPATH.'%(genericUtils_path)
sys.path.append(genericUtils_path)

ROOT.gROOT.SetMacroPath('%s:%s/share'%(ROOT.gROOT.GetMacroPath(),the_path))
ROOT.gROOT.LoadMacro('CutsMacro.h')

ROOT.gROOT.SetBatch(True)

import python.PlotFunctions as plotfunc
import python.TAxisFunctions as taxisfunc
import python.PyAnalysisPlotting as anaplot
from python.ShowerShapeEvolutionPlotter import ShowerShapeEvolutionPlot
import PhotonIDHelpers as idhelpers

#-----------------------------------------------
def main(options,args) :

    plotfunc.SetupStyle()
    options.names = options.names.split(',')

    confs = dict()
    for confstr in ['tight','loose','menu3','menu4'] :
        print 'Using %s for %s'%(getattr(options,confstr),confstr)
        confs[confstr] = ROOT.TEnv(getattr(options,confstr))

    # Assume for now that Strips and non-strips binning is the same
    eta_bins = [0,0.6,0.8,1.15,1.37,1.52,1.81,2.01,2.37,2.47]
    print eta_bins

    # Assume for now that Strips and non-strips binning is the same
    et_bins = [10000,15000,20000,25000,30000,35000,40000,45000,50000,60000,80000,100000]
    print et_bins

    # Make dummy histograms, because the "EvaluatePhotonID" function needs it
    args = len(et_bins)-1,array('f',list(a/1000. for a in et_bins)),len(eta_bins)-1,array('f',eta_bins)
    numerator_tmp   = ROOT.TH2F('DUMMY_num','DUMMY',*args)
    denominator_tmp = ROOT.TH2F('DUMMY_den','DUMMY',*args)

    if options.makehists :

        outputfile = ROOT.TFile('output.root','RECREATE');

        for status in ['Converted','NonConverted'] :

            print 'Making plots for %s'%(status)

            trees = dict()

            files_rz ,trees_rz ,keys_rz  = anaplot.GetTreesFromFiles(options.radzsignal        ,treename=options.radztreename        )
            files_rzd,trees_rzd,keys_rzd = anaplot.GetTreesFromFiles(options.radzdata          ,treename=options.radztreename        )
            files_sp ,trees_sp ,keys_sp  = anaplot.GetTreesFromFiles(options.singlephotonsignal,treename=options.singlephotontreename)
            files_spd,trees_spd,keys_spd = anaplot.GetTreesFromFiles(options.singlephotondata  ,treename=options.singlephotontreename)
            files_jf ,trees_jf ,keys_jf  = anaplot.GetTreesFromFiles(options.jetfiltered       ,treename=options.singlephotontreename)

            trees['radz']      = trees_rz [keys_rz [0]] if keys_rz  else None
            trees['radz_data'] = trees_rzd[keys_rzd[0]] if keys_rzd else None
            trees['incl']      = trees_sp [keys_sp [0]] if keys_sp  else None
            trees['incl_data'] = trees_spd[keys_spd[0]] if keys_spd else None
            trees['jf']        = trees_jf [keys_jf [0]] if keys_jf  else None

            variables = dict()
            variables['radz']      = ROOT.photonVariables('RadZ_%s'             %(status),numerator_tmp.GetNbinsX(),numerator_tmp.GetNbinsY(),status == 'Converted')
            variables['radz_data'] = ROOT.photonVariables('RadZ_data_%s'        %(status),numerator_tmp.GetNbinsX(),numerator_tmp.GetNbinsY(),status == 'Converted')
            variables['incl']      = ROOT.photonVariables('SinglePhoton_%s'     %(status),numerator_tmp.GetNbinsX(),numerator_tmp.GetNbinsY(),status == 'Converted')
            variables['incl_data'] = ROOT.photonVariables('SinglePhoton_data_%s'%(status),numerator_tmp.GetNbinsX(),numerator_tmp.GetNbinsY(),status == 'Converted')
            variables['jf']        = ROOT.photonVariables('Jetfiltered_%s'      %(status),numerator_tmp.GetNbinsX(),numerator_tmp.GetNbinsY(),status == 'Converted')

            doSample = dict()
            doSample['radz']      = options.radzsignal
            doSample['radz_data'] = options.radzdata
            doSample['incl']      = options.singlephotonsignal
            doSample['incl_data'] = options.singlephotondata
            doSample['jf']        = options.jetfiltered

            colors = dict()
            colors['radz']      = ROOT.kAzure-2
            colors['radz_data'] = ROOT.kBlack
            colors['incl']      = ROOT.kRed+1
            colors['incl_data'] = ROOT.kGray+2 if doSample['radz_data'] else ROOT.kBlack
            colors['jf']        = ROOT.kGreen+1

            titles = dict()
            titles['radz']      = 'Z#rightarrow^{}ll#gamma MC (fudged)'
            titles['radz_data'] = 'Z#rightarrow^{}ll#gamma data'
            titles['incl']      = 'Inclusive photons MC (fudged)'
            titles['incl_data'] = 'Inclusive photons data'
            titles['jf']        = 'Jet-filtered MC (fudged)'

            # general setup of tight ID menu
            etbins_tight  = idhelpers.GetCutValuesFromConf(confs['tight'],'CutBinEnergy',status)
            etabins_tight =  idhelpers.GetCutValuesFromConf(confs['tight'],'CutBinEta',status)
            n_et_tight  = 1 + len(etbins_tight) # see note above
            n_eta_tight = len(etabins_tight)
            tight_id = ROOT.photonID(n_et_tight,n_eta_tight)
            tight_id.Set_EtaBinThresholds(array('f',etabins_tight))
            if etbins_tight :
                tight_id.Set_EtBinThresholds(array('f',etbins_tight))

            for var in list(variables['radz'].variables) :
                cut_values = idhelpers.GetCutValuesFromConf(confs['tight'],var,status)
                cut_values = array('f',cut_values)
                getattr(tight_id,'Set_%s'%(var))(cut_values)

            # Evaluate photon ID for each of the samples
            for sample in trees.keys() :
                if not doSample[sample] :
                    continue
                print 'Evaluating ID for %s'%(sample)

                doTruthMatchPhoton = {'radz'     :False,
                                      'radz_data':False,
                                      'incl'     :True,
                                      'incl_data':False,
                                      'jf'       :False,
                                      }.get(sample)
                doTruthMatchFake = {'radz'     :False,
                                    'radz_data':False,
                                    'incl'     :False,
                                    'incl_data':False,
                                    'jf'       :True,
                                    }.get(sample)

                ROOT.EvaluatePhotonID(trees[sample],
                                      tight_id,
                                      status == 'Converted',
                                      denominator_tmp,
                                      numerator_tmp,
                                      variables[sample],
                                      doTruthMatchPhoton,
                                      doTruthMatchFake,
                                      options.FixedCutLoose)

            # Process the cuts and the histograms
            for var in list(variables['radz'].variables) :

                # Reset PlotFunctions tobject_collector
                del plotfunc.tobject_collector[:]

                for eta in range(len(eta_bins)-1) :

                    if eta_bins[eta] in [1.37,2.37] :
                        continue

                    hists = dict()
                    for sample in trees.keys() :
                        hists[sample] = []

                    cuts_graphs = dict()

                    names = ['Tight','Loose','Menu3','Menu4']
                    for ii in range(len(options.names)) :
                        names[ii] = options.names[ii]

                    for et in range(len(et_bins)-1) :

                        for sample in trees.keys() :
                            if not doSample[sample] :
                                continue
                            hist = getattr(variables[sample],'h_%s'%(var))[et*variables[sample].neta + eta]
                            hist.SetTitle(titles[sample])
                            hist.SetLineColor(colors[sample])
                            hist.SetLineWidth(2)
                            if 'data' in sample :
                                hist.SetMarkerSize(0.5)
                                hist.SetOption('pE')
                            else :
                                hist.SetOption('hist')

                            outputfile.cd()
                            hist.Write()

            for i in files_rz.keys() :
                files_rz[i].Close()
            for i in files_sp.keys() :
                files_sp[i].Close()

            ## end variable block
        ## end status block

        outputfile.Close()


    
    if options.makefits :

        outputfile = ROOT.TFile('output.root','READ');
        hist = outputfile.Get('CutHadLeakage_SinglePhoton_Converted_3_0')

        ROOT.gROOT.SetBatch(False)
#         hist.Draw()

#         w = ROOT.RooWorkspace("w",ROOT.kTRUE)
#         func = w.factory("Gaussian::gauss(variable[-0.05,0.05],mean[5.28,5.2,5.3],width[0.0027,0.001,1])")

        w = ROOT.RooWorkspace("w",ROOT.kTRUE)
        name_variable = 'rhad'
        min = hist.GetBinLowEdge(0)
        max = hist.GetBinLowEdge(hist.GetNbinsX()+1)
        #variable = ROOT.RooRealVar(name_variable,name_variable,min,max)
        variable = w.factory('%s[%f,%f]'%(name_variable,min,max))

        name_function = 'decay'

        if name_function == 'gaussian' :
            sigmean  = ROOT.RooRealVar("sigmean" ,"mass"  ,0,-0.05,0.05)
            sigwidth = ROOT.RooRealVar("sigwidth","width" ,0.01,0.0001,1.)
            func = ROOT.RooGaussian("func","func",variable,sigmean,sigwidth)
    
        elif name_function == 'double gaussian' :
            w.factory("Gaussian::g1(%s,mean_g1[0,-0.05,0.05],sigma_g1[0.01,0,0.1])"%(name_variable))
            w.factory("Gaussian::g2(%s,mean_g2[0.03,-0.05,0.05],sigma_g2[0.01,0,0.1])"%(name_variable))
            func = w.factory("SUM::model(frac[0.5,0,1]*g1,g2)") ;
        
        elif name_function == 'triple gaussian' :
            w.factory("Gaussian::g1(%s,mean_g1[0,-0.05,0.05],sigma_g1[0.01,0,0.1])"%(name_variable))
            w.factory("Gaussian::g2(%s,mean_g2[0.03,-0.05,0.05],sigma_g2[0.01,0,0.1])"%(name_variable))
            w.factory("Gaussian::g3(%s,mean_g3[0.03,-0.05,0.05],sigma_g3[0.01,0,0.1])"%(name_variable))
            func = w.factory("SUM::model(frac_g1[0.5,0,1]*g1,frac_g2[0.5,0,1]*g2,g3)") ;
        
        elif name_function == 'quadruple gaussian' :
            w.factory("Gaussian::g1(%s,mean_g1[0,-0.05,0.05],sigma_g1[0.01,0,0.1])"%(name_variable))
            w.factory("Gaussian::g2(%s,mean_g2[0,-0.05,0.05],sigma_g2[0.01,0,0.1])"%(name_variable))
            w.factory("Gaussian::g3(%s,mean_g3[0,-0.05,0.05],sigma_g3[0.01,0,0.1])"%(name_variable))
            w.factory("Gaussian::g4(%s,mean_g4[0,-0.05,0.05],sigma_g4[0.01,0,0.1])"%(name_variable))
            func = w.factory("SUM::model(frac_g1[0.5,0,1]*g1,frac_g2[0.5,0,1]*g2,frac_g3[0.5,0,1]*g3,g4)") ;
        
        elif name_function == 'breitwigner' :
            func = w.factory("BreitWigner::g1(%s,mean_g1[0,-0.05,0.05],sigma_g1[0.01,0,0.1])"%(name_variable))

        elif name_function == 'bw+gaus' :
            w.factory("BreitWigner::b1(%s,mean_b1[0,-0.05,0.05],sigma_b1[0.01,0,0.1])"%(name_variable))
            w.factory("Gaussian::g1(%s,mean_g1[0,-0.05,0.05],sigma_g1[0.01,0,0.1])"%(name_variable))
            func = w.factory("SUM::model(frac[0.5,0,1]*g1,b1)") ;

        else :
            w.factory("BreitWigner::b1(%s,mean_b1[0,-0.05,0.05],sigma_b1[0.01,0,0.1])"%(name_variable))
            w.factory("Gaussian::g1(%s,mean_g1[0,-0.05,0.05],sigma_g1[0.01,0,0.1])"%(name_variable))
            w.factory("Gaussian::g2(%s,mean_g2[0,-0.05,0.05],sigma_g2[0.01,0,0.1])"%(name_variable))
            func = w.factory("SUM::model(frac[0.5,0,1]*g1,frac_g2[0.5,0,1]*g2,b1)") ;

        data = ROOT.RooDataHist("data","data",ROOT.RooArgList(variable),hist);

        func.fitTo(data)
        frame = variable.frame()
        data.plotOn(frame)
        func.plotOn(frame)
        
#         frame.Draw()
#         raw_input('pause')

        roofit_hist = frame.getHist('h_data')
        curve = frame.getCurve()
        pull = roofit_hist.makePullHist(curve)
        
        can = plotfunc.RatioCanvas('can','can')
        plotfunc.AddHistogram(can,roofit_hist)
#         plotfunc.AddHistogram(can,curve,'l')
#        plotfunc.AddRatio(can,curve,roofit_hist)
        plotfunc.AddRatioManual(can,curve,pull,drawopt1='l',drawopt2='pE1')
        plotfunc.FormatCanvasAxes(can)

        raw_input('pause')

        outputfile.Close()

    return
    
#-----------------------------------------------
if __name__ == '__main__':
    from optparse import OptionParser
    p = OptionParser()
    p.add_option('--tight',type = 'string', default = '', dest = 'tight',help = 'Tight Menu')
    p.add_option('--loose' ,type = 'string', default = '', dest = 'loose' ,help = 'Loose Menu' )
    p.add_option('--menu3' ,type = 'string', default = '', dest = 'menu3' ,help = 'Another Menu' )
    p.add_option('--menu4' ,type = 'string', default = '', dest = 'menu4' ,help = 'Another Menu' )

    p.add_option('--names',type = 'string', default = '', dest = 'names',help = 'Menu names, comma-separated (tight,loose,menu3,menu4')

    p.add_option('--radzsignal'  ,type = 'string', default = '', dest = 'radzsignal' ,help = 'Radiative-Z Signal file' )
    p.add_option('--radzdata'    ,type = 'string', default = '', dest = 'radzdata'   ,help = 'Radiative-Z Data file'   )
    p.add_option('--radztreename',type = 'string', default = 'output', dest = 'radztreename' ,help = 'Radiative-Z treename' )

    p.add_option('--singlephotonsignal'  ,type = 'string', default = '', dest = 'singlephotonsignal' ,help = 'Single photon Signal file' )
    p.add_option('--singlephotondata'    ,type = 'string', default = '', dest = 'singlephotondata'   ,help = 'Single photon Data file'   )
    p.add_option('--jetfiltered'         ,type = 'string', default = '', dest = 'jetfiltered'        ,help = 'Jet-filtered MC file'      )
    p.add_option('--singlephotontreename',type = 'string', default = 'SinglePhoton', dest = 'singlephotontreename' ,help = 'Single photon treename' )

    p.add_option('--FixedCutLoose',action='store_true',default=False,dest='FixedCutLoose',help='apply FixedCutLoose preselection')

    p.add_option('--makehists',action='store_true',default=False,dest='makehists',help='Make the histograms')
    p.add_option('--makefits' ,action='store_true',default=False,dest='makefits' ,help='Do the fits')

    p.add_option('--outdir',type='string',default='.',dest='outdir',help='output directory')

    options,args = p.parse_args()
    
    main(options,args)
