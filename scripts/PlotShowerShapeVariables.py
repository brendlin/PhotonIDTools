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

            can_barrel = ROOT.TCanvas('can_%s_%s_barrel'%(var,status),'blah',int(1464),959)
            can_endcap = ROOT.TCanvas('can_%s_%s_endcap'%(var,status),'blah',int(1464),959)
            pads = []

            for eta in range(len(eta_bins)-1) :

                if eta_bins[eta] in [1.37,2.37] :
                    continue

                # Bin labeling
                template = '%2.0f^{ }<^{ }p_{T}^{ }<^{ }%2.0f GeV'
                labels = list(template%(et_bins[a]/1000,et_bins[a+1]/1000) for a in range(len(et_bins)-1))

                hists = dict()
                for sample in trees.keys() :
                    hists[sample] = []

                def GetCutsGraph(gr_title) :
                    tmp_gr = ROOT.TGraph(2*len(et_bins)-2,
                                         array('f',[0]*(2*len(et_bins)-2)),
                                         array('f',[0]*(2*len(et_bins)-2)))
                    tmp_gr.SetTitle(gr_title)
                    tmp_gr.SetName('%s_%s_%s_%d'%(gr_title.replace(' ','_'),var,status,eta))
                    tmp_gr.SetLineWidth(1)
                    return tmp_gr

                cuts_graphs = dict()

                names = ['Tight','Loose','Menu3','Menu4']
                for ii in range(len(options.names)) :
                    names[ii] = options.names[ii]

                # Tight cuts graph
                cuts_graphs['tight'] = GetCutsGraph(names[0])
                cuts_graphs['tight'].SetLineColor(ROOT.kBlue+1)

                # Loose cuts graph
                cuts_graphs['loose'] = GetCutsGraph(names[1])
                cuts_graphs['loose'].SetLineColor(ROOT.kOrange+5)
                cuts_graphs['loose'].SetLineWidth(2)

                # menu 3 graph
                cuts_graphs['menu3'] = GetCutsGraph(names[2])
                cuts_graphs['menu3'].SetLineColor(ROOT.kGray+2)
                cuts_graphs['menu3'].SetLineStyle(2)

                # menu 4 graph
                cuts_graphs['menu4'] = GetCutsGraph(names[3])
                cuts_graphs['menu4'].SetLineColor(ROOT.kGreen+1)
                cuts_graphs['menu4'].SetLineStyle(7)

                for et in range(len(et_bins)-1) :

                    for sample in trees.keys() :
                        if not doSample[sample] :
                            continue
                        hists[sample].append(getattr(variables[sample],'h_%s'%(var))[et*variables[sample].neta + eta])
                        hists[sample][-1].SetTitle(titles[sample])
                        hists[sample][-1].SetLineColor(colors[sample])
                        hists[sample][-1].SetLineWidth(2)
                        if 'data' in sample :
                            hists[sample][-1].SetMarkerSize(0.5)
                            hists[sample][-1].SetOption('pE')
                        else :
                            hists[sample][-1].SetOption('hist')

                    # Get the cut values
                    for id in ['tight','loose','menu3','menu4'] :

                        # Get the right cut threshold
                        # Note that the thresholds must match between the plotting and the cut menu!
                        et_cutvalues = et
                        etbins_tight = idhelpers.GetCutValuesFromConf(confs[id],'CutBinEnergy',status)

                        if etbins_tight :
                            etbins_tight  = [0] + list(int(a) for a in etbins_tight)
                            if et_bins[et] > etbins_tight[-1] :
                                et_cutvalues = len(etbins_tight)-1
                            elif et_bins[et] < etbins_tight[1] :
                                # bins between 0 and the first bin
                                et_cutvalues = 0
                            else :
                                et_cutvalues = etbins_tight.index(int(et_bins[et]))

                        cut_value = idhelpers.GetCutValueFromConf(confs[id],var,status,et_cutvalues,eta)
                        if cut_value == None :
                            continue
                        cuts_graphs[id].SetPoint(et*2  ,cut_value,et  )
                        cuts_graphs[id].SetPoint(et*2+1,cut_value,et+1)

                # some plot text
                text_lines = [plotfunc.GetAtlasInternalText()]
                text_lines += [plotfunc.GetSqrtsText(13)+', 13 fb^{#minus1}']
                text_lines += ['%2.2f < |#eta| < %2.2f'%(eta_bins[eta],eta_bins[eta+1])]

                # Individual canvas
                if False :
                    can = ROOT.TCanvas('can_%s_%s_%s'%(var,status,eta),'blah',520,700)
                    ShowerShapeEvolutionPlot(can,labels,hists['radz'],hists['incl'],hists['jf'],hists['radz_data'],hists['incl_data'])
                    plotfunc.AddHistogram(can,cuts_graphs['loose'],drawopt='l')
                    plotfunc.AddHistogram(can,cuts_graphs['tight'],drawopt='l')
                    if options.menu3 :
                        plotfunc.AddHistogram(can,cuts_graphs['menu3'],drawopt='l')
                    if options.menu4 :
                        plotfunc.AddHistogram(can,cuts_graphs['menu4'],drawopt='l')
                    plotfunc.SetLeftMargin(can,0.33)
                    plotfunc.MakeLegend(can,0.7,0.88,0.9,0.99,totalentries=3)

                    plotfunc.DrawText(can,text_lines,0.33,0.88,0.5,0.99,totalentries=3)

                    can.SetGridx()
                    taxisfunc.SetNdivisions(can,5,5,0)
                    can.Update()

                    if not os.path.exists(status) :
                        os.makedirs(status)
                    can.Print('%s/%s_%s_%s.pdf'%(status,var,status,eta))

                # Full barrel and endcap canvases
                if True :
                    can_barrel.cd() if (eta < 5) else can_endcap.cd()
                    eta1 = eta%5
                    pads.append(ROOT.TPad('pad_%s_%s_%s'%(var,status,eta),'blah',(172*(eta1>0) + eta1*348)/1564.,0,(172 + 348*(eta1+1))/1564.,1))
                    ShowerShapeEvolutionPlot(pads[-1],labels,hists['radz'],hists['incl'],hists['jf'],hists['radz_data'],hists['incl_data'])
                    plotfunc.AddHistogram(pads[-1],cuts_graphs['loose'],drawopt='l')
                    plotfunc.AddHistogram(pads[-1],cuts_graphs['tight'],drawopt='l')
                    if options.menu3 :
                        plotfunc.AddHistogram(pads[-1],cuts_graphs['menu3'],drawopt='l')
                    if options.menu4 :
                        plotfunc.AddHistogram(pads[-1],cuts_graphs['menu4'],drawopt='l')

                    pads[-1].SetBottomMargin(0.06)
                    if not eta1 :
                        plotfunc.SetLeftMargin(pads[-1],0.33)
                    if eta1 :
                        plotfunc.SetLeftMargin(pads[-1],0.01)
                    if eta1 < 3 :
                        plotfunc.SetRightMargin(pads[-1],0.01)
                    can_barrel.cd() if (eta < 5) else can_endcap.cd()
                    pads[-1].Draw()
                    taxisfunc.SetNdivisions(pads[-1],5,5,0)
                    for i in pads[-1].GetListOfPrimitives() :
                        if hasattr(i,'SetLineWidth') and issubclass(type(i),ROOT.TH1) :
                            i.SetLineWidth(1)
                        if hasattr(i,'GetYaxis') :
                            i.GetYaxis().SetLabelSize(19)
                            i.GetXaxis().SetLabelSize(19)
                        if issubclass(type(i),ROOT.TH1) and not eta1 :
                            i.SetTitle('remove')
                        if issubclass(type(i),ROOT.TGraph) and eta1 == 1 :
                            i.SetTitle('remove')

                    if not eta1 :
                        plotfunc.MakeLegend(pads[-1],0.7,0.88,0.9,0.99,totalentries=3)
                        pads[-1].GetPrimitive('legend').SetMargin(0.4)
                        plotfunc.DrawText(pads[-1],text_lines,0.33,0.88,0.5,0.99,totalentries=3)
                    else :
                        plotfunc.DrawText(pads[-1],['','',text_lines[-1]],0.01,0.88,0.5,0.99,totalentries=1)
                        if pads[-1].GetPrimitive('legend') :
                            pads[-1].GetPrimitive('legend').Delete()

                    if eta1 == 1 :
                        plotfunc.MakeLegend(pads[-1],0.01,0.93,0.5,0.99,totalentries=2)

                ## end barrel / EC block
            ## end eta block

            if not os.path.exists('%s/%s'%(options.outdir,status)) :
                os.makedirs('%s/%s'%(options.outdir,status))

            can_barrel.Print('%s/%s/barrel_%s_%s.pdf'%(options.outdir,status,var,status))
            can_endcap.Print('%s/%s/endcap_%s_%s.pdf'%(options.outdir,status,var,status))

            # Memory management
            for pad in list(can_barrel.GetListOfPrimitives()) :
                if not issubclass(type(pad),ROOT.TPad) :
                    continue
                for hist in list(pad.GetListOfPrimitives()) :
                    if not issubclass(type(hist),ROOT.TH1) :
                        continue
                    hist.Delete()
            ROOT.gROOT.ProcessLine('delete %s'%(can_barrel.GetName()))
            ROOT.gROOT.ProcessLine('delete %s'%(can_endcap.GetName()))

        for i in files_rz.keys() :
            files_rz[i].Close()
        for i in files_sp.keys() :
            files_sp[i].Close()

        ## end variable block
    ## end status block
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

    p.add_option('--outdir',type='string',default='.',dest='outdir',help='output directory')

    options,args = p.parse_args()
    
    main(options,args)
