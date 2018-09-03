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

##################################################################################
def getNdof(thing) :
    argSet = thing.getVariables()
    iter = argSet.createIterator()
    var = iter.Next()
    nNotConstant = 0
    while var :
        if not var.isConstant() :
            nNotConstant += 1
        var = iter.Next()
    return nNotConstant

##################################################################################
def printArgSetArgs(argSet,nest='') :
    iter = argSet.createIterator()
    var = iter.Next()
    while var :
        print nest+' ',var.GetName(),type(var),var.getVal(),'isConstant:',var.isConstant(),
        if hasattr(var,'getMin') :
            print 'min',var.getMin(),'max',var.getMax()
        else :
            print
        var = iter.Next()
    return

##################################################################################
def PrintRooThing(thing,nest='') :
    print nest+thing.GetName(),type(thing)

    if hasattr(thing,'getComponents') : # a RooAbsArg function
        print nest+'Components:'
        printArgSetArgs(thing.getComponents(),nest=nest)
        iter = thing.getComponents().createIterator()
        var = iter.Next()
        while var :
            if var == thing :
                var = iter.Next()
                continue
            PrintRooThing(var,nest=nest+'   ')
            var = iter.Next()
    if hasattr(thing,'getVariables') : # a RooAbsArg function
        print nest+'Variables:'
        printArgSetArgs(thing.getVariables(),nest=nest)

    return


#-----------------------------------------------
def GetGaussianWithShiftAndWidthOptions(name,workspace,min,max,transformname = '') :

    average = (min+max)/2.
    step_size = (max-min)/10.
    if not workspace.var('x') :
        workspace.factory('x[%f,%f,%f]'%(average,min,max))
    if not workspace.var('shift_%s'%(transformname)) :
        workspace.factory('shift_%s[%f,%f,%f]'%(transformname,0,-step_size,step_size))
    if not workspace.var('width_%s'%(transformname)) :
        workspace.factory('width_%s[%f,%f,%f]'%(transformname,1,0,2))

    muname = 'mu_%s'%(name)
    sigmaname = 'sigma_%s'%(name)

    mu = workspace.factory('%s[%f,%f,%f]'%(muname,average,min,max))
    sigma = workspace.factory('%s[%f,%f,%f]'%(sigmaname,step_size,step_size/5.,5*step_size))

    x_shifted = 'x-shift_%s'%(transformname)
    x_shifted_widened = '((x-shift_%s)/width_%s)'%(transformname,transformname)
    x_widened = '(x/width_%s)'%(transformname)

    args = (x_shifted_widened,muname,x_shifted_widened,muname,sigmaname,sigmaname)
    gauss_expr = 'exp(-0.5*(%s-%s)*(%s-%s)/(%s*%s))'%(args)
    ws_expression = "EXPR::%s('%s',x,shift_%s,width_%s,%s,%s)"%(name,gauss_expr,transformname,transformname,muname,sigmaname)
    gauss = workspace.factory(ws_expression)

    workspace.var('shift_%s'%(transformname)).setConstant(ROOT.kTRUE)
    workspace.var('width_%s'%(transformname)).setConstant(ROOT.kTRUE)

    return gauss

#-----------------------------------------------
def GetNGaussian(nGaus,name,workspace,min,max) :
    transformname = 'Vg'
    rndm = ROOT.TRandom3()
    average = (min+max)/2.
    for i in range(nGaus) :
        GetGaussianWithShiftAndWidthOptions('%s_%s%d'%(name,transformname,i),workspace,min,max,transformname=transformname)
        muval = average*(0.8 + rndm.Rndm()*0.4)
        sigval = workspace.var('sigma_%s_%s%d'%(name,transformname,i)).getVal() * (1+rndm.Rndm()*1.0)
        #print 'muval:',muval,'sigval:',sigval
        workspace.var('mu_%s_%s%d'%(name,transformname,i)).setVal(muval)
        workspace.var('sigma_%s_%s%d'%(name,transformname,i)).setVal(sigval)

    fac_string = 'SUM::%s ( %s_%s0'%(name,name,transformname)
    for i in range(1,nGaus) :
        start = i/float(nGaus)
        fac_string += ', frac%d_%s_%s[%f,0.001,1]*'%(i,transformname,name,start)
        fac_string += '%s_%s%d'%(name,transformname,i)
    fac_string += ' )'
    print fac_string
    Ngauss = workspace.factory(fac_string)

    return Ngauss

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


    #
    # Step 2: fit signal and background.
    #
    if options.makefits :

        outputfile = ROOT.TFile('output.root','READ');
        hist = outputfile.Get('fracm_SinglePhoton_Converted_3_0')
        hbkg = outputfile.Get('fracm_Jetfiltered_Converted_3_0')
        hdata = outputfile.Get('fracm_SinglePhoton_data_Converted_3_0')

        ROOT.gROOT.SetBatch(False)

        w = ROOT.RooWorkspace("w",ROOT.kTRUE)
        name_variable = 'x'
        min = hist.GetBinLowEdge(5)
        max = hist.GetBinLowEdge(hist.GetNbinsX()+1)
        variable = w.factory('%s[%f,%f]'%(name_variable,min,max))

        fsig = GetNGaussian(4,'sig',w,min,max)
        fbkg = GetNGaussian(4,'bkg',w,min,max)

        mc   = ROOT.RooDataHist("mc"  ,"mc"  ,ROOT.RooArgList(variable),hist)
        bkg  = ROOT.RooDataHist('bkg' ,'bkg' ,ROOT.RooArgList(variable),hbkg)
        data = ROOT.RooDataHist('data','data',ROOT.RooArgList(variable),hdata)

        frame = variable.frame()
        mc.plotOn(frame)
        mc_hist = frame.getHist('h_mc')

        can = plotfunc.RatioCanvas('can','can')
        plotfunc.AddHistogram(can,mc_hist)

        doPull = False

        fsig.fitTo(mc)
        fsig.plotOn(frame,ROOT.RooFit.Name(fsig.GetName()))
        curve = frame.getCurve(fsig.GetName())
        pull = mc_hist.makePullHist(curve) if doPull else mc_hist.makeResidHist(curve)
        PrintRooThing(fsig)
        ndof = getNdof(fsig) - 1 # subtract one for "x"
        chisquare = frame.chiSquare(1+ndof) # n-1 bins
        bins = hist.GetNbinsX()
        pvalue_chi2 = ROOT.TMath.Prob(chisquare*(bins-1-ndof),bins-1-ndof)
        #print 'ndof:',ndof,'chisquare:',chisquare,'pvalue:',pvalue_chi2
        plotfunc.AddRatioManual(can,curve,pull,drawopt1='l',drawopt2='pE1')

        bframe = variable.frame(ROOT.RooFit.Title("bkg distribution"))
        bkg.plotOn(bframe)
        bkg_hist = bframe.getHist('h_bkg')

        bcan = plotfunc.RatioCanvas('bcan','bcan')
        plotfunc.AddHistogram(bcan,bkg_hist)

        fbkg.fitTo(bkg)
        fbkg.plotOn(bframe,ROOT.RooFit.Name(fbkg.GetName()))
        curve = bframe.getCurve(fbkg.GetName())
        pull = bkg_hist.makePullHist(curve) if doPull else bkg_hist.makeResidHist(curve)
        plotfunc.AddRatioManual(bcan,curve,pull,drawopt1='l',drawopt2='pE1')

        dframe = variable.frame(ROOT.RooFit.Title("data distribution"))
        data.plotOn(dframe)
        data_hist = dframe.getHist('h_data')

        dcan = plotfunc.RatioCanvas('dcan','dcan')
        plotfunc.AddHistogram(dcan,data_hist)

        # fix model
        for f in [fsig,fbkg] :
            argSet = f.getVariables()
            iter = argSet.createIterator()
            var = iter.Next()
            while var :
                name = var.GetName()
                if name != 'x' and 'shift' not in name :
                    var.setConstant(ROOT.kTRUE)
                if 'shift' in name :
                    var.setConstant(ROOT.kFALSE)
                #print var.GetName(),var.isConstant()
                var = iter.Next()

        model = w.factory('SUM::model ( frac_model[0.5,0,1]*%s, %s)'%(fsig.GetName(),fbkg.GetName()))
        model.fitTo(data)
        model.plotOn(dframe,ROOT.RooFit.Name(model.GetName()))
        model.plotOn(dframe,ROOT.RooFit.Components(fbkg.GetName()),ROOT.RooFit.Name('model_bkg_component'),ROOT.RooFit.LineStyle(ROOT.kDashed))
        model.plotOn(dframe,ROOT.RooFit.Components(fsig.GetName()),ROOT.RooFit.Name('model_sig_component'))

        curve = dframe.getCurve(model.GetName())
        pull = data_hist.makePullHist(curve) if doPull else data_hist.makeResidHist(curve)
        plotfunc.AddRatioManual(dcan,curve,pull,drawopt1='l',drawopt2='pE1')
        plotfunc.AddHistogram(dcan,dframe.getCurve('model_bkg_component'),'l')
        plotfunc.AddHistogram(dcan,dframe.getCurve('model_sig_component'),'l')

        # print model
        for f in [model] :
            argSet = f.getVariables()
            iter = argSet.createIterator()
            var = iter.Next()
            while var :
                print var.GetName(),var.isConstant(),var.getVal()
                var = iter.Next()

        for c in [can,bcan,dcan] :
            plotfunc.SetColors(c)
            plotfunc.FormatCanvasAxes(c)

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
