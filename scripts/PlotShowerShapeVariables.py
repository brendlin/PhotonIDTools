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

    for status in ['Converted','NonConverted'] :

        print 'Making plots for %s'%(status)

        files_rz,trees_rz,keys_rz = anaplot.GetTreesFromFiles(options.radzsignal,treename=options.radztreename)
        files_sp,trees_sp,keys_sp = anaplot.GetTreesFromFiles(options.singlephotonsignal,treename=options.singlephotontreename)

        # Assume for now that Strips and non-strips binning is the same
        eta_bins = confs['tight'].GetValue('CutBinEta_photons%s'%(status),'').split(';')
        eta_bins = list(float(a.rstrip().lstrip()) for a in eta_bins)
        eta_bins = [0] + eta_bins
        print eta_bins

        # Assume for now that Strips and non-strips binning is the same
        et_bins = confs['tight'].GetValue('CutBinEnergy_photons%s'%(status),'').split(';')
        et_bins = list(float(a.rstrip().lstrip()) for a in et_bins)
        et_bins = [10000] + et_bins
        tmp = []
        for i in et_bins :
            if i < 100000.1 :
                tmp.append(i)
        et_bins = tmp
        print et_bins

        args = len(et_bins)-1,array('d',list(a/1000. for a in et_bins)),len(eta_bins)-1,array('d',eta_bins)
        numerator_rz   = ROOT.TH2F('numerator_%s_%s'  %(status,'radz'),'Z#rightarrow^{}ll#gamma MC (fudged)',*args)
        denominator_rz = ROOT.TH2F('denominator_%s_%s'%(status,'radz'),'Z#rightarrow^{}ll#gamma MC (fudged)',*args)

        numerator_sp   = ROOT.TH2F('numerator_%s_%s'  %(status,'incl'),'Inclusive #gamma MC (fudged)'       ,*args)
        denominator_sp = ROOT.TH2F('denominator_%s_%s'%(status,'incl'),'Inclusive #gamma MC (fudged)'       ,*args)

        variables_rz = ROOT.photonVariables('RadZ_%s'        %(status),numerator_rz.GetNbinsX(),numerator_rz.GetNbinsY(),status == 'Converted')
        variables_sp = ROOT.photonVariables('SinglePhoton_%s'%(status),numerator_sp.GetNbinsX(),numerator_sp.GetNbinsY(),status == 'Converted')

        # general setup of tight ID menu
        etbins_tight  = idhelpers.GetCutValuesFromConf(confs['tight'],'CutBinEnergy',status)
        etabins_tight =  idhelpers.GetCutValuesFromConf(confs['tight'],'CutBinEta',status)
        n_et_tight  = 1 + len(etbins_tight) # see note above
        n_eta_tight = len(etabins_tight)
        tight_id = ROOT.photonID(n_et_tight,n_eta_tight)
        tight_id.Set_EtaBinThresholds(array('d',etabins_tight))
        tight_id.Set_EtBinThresholds(array('d',etbins_tight))

        for var in list(variables_rz.variables) :
            cut_values = idhelpers.GetCutValuesFromConf(confs['tight'],var,status)
            cut_values = array('d',cut_values)
            getattr(tight_id,'Set_%s'%(var))(cut_values)

        # New Radiative-Z efficiency method:
        if options.radzsignal :
            print 'Evaluating Radiative-Z signal ID...'
            ROOT.EvaluatePhotonID_InclusivePhoton(trees_rz[keys_rz[0]],
                                                  tight_id,
                                                  status == 'Converted',
                                                  denominator_rz,
                                                  numerator_rz,
                                                  variables_rz)

        # New Single Photon efficiency method:
        if options.singlephotonsignal :
            print 'Evaluating Single-photon signal ID...'
            ROOT.EvaluatePhotonID_InclusivePhoton(trees_sp[keys_sp[0]],
                                                  tight_id,
                                                  status == 'Converted',
                                                  denominator_sp,
                                                  numerator_sp,
                                                  variables_sp)


        # Process the cuts and the histograms
        for var in list(variables_rz.variables) :

            # Reset PlotFunctions tobject_collector
            del plotfunc.tobject_collector[:]

            can_barrel = ROOT.TCanvas('can_%s_%s_barrel'%(var,status),'blah',int(1464),700)
            can_endcap = ROOT.TCanvas('can_%s_%s_endcap'%(var,status),'blah',int(1464),700)
            pads = []

            for eta in range(len(eta_bins)-1) :

                if eta_bins[eta] in [1.37,2.37] :
                    continue

                # Bin labeling
                template = '%2.0f^{ }<^{ }p_{T}^{ }<^{ }%2.0f GeV'
                labels = list(template%(et_bins[a]/1000,et_bins[a+1]/1000) for a in range(len(et_bins)-1))

                radz_hists = []
                singlephoton_hists = []

                def GetCutsGraph(gr_title) :
                    tmp_gr = ROOT.TGraph(2*len(et_bins)-2,
                                         array('d',[0]*(2*len(et_bins)-2)),
                                         array('d',[0]*(2*len(et_bins)-2)))
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

                    if options.radzsignal :
                        radz_hists.append(getattr(variables_rz,'h_%s'%(var))[et*variables_rz.neta + eta])
                        radz_hists[-1].SetTitle('Z#rightarrow^{}ll#gamma MC (fudged)')
                        radz_hists[-1].SetLineColor(ROOT.kBlack)
                        radz_hists[-1].SetLineWidth(2)

                    if options.singlephotonsignal :
                        singlephoton_hists.append(getattr(variables_sp,'h_%s'%(var))[et*variables_sp.neta + eta])
                        singlephoton_hists[-1].SetTitle('Inclusive photons (fudged)')
                        singlephoton_hists[-1].SetLineColor(ROOT.kRed+1)
                        singlephoton_hists[-1].SetLineWidth(2)

                    # Get the cut values
                    for id in ['tight','loose','menu3','menu4'] :
                        cut_value = idhelpers.GetCutValueFromConf(confs[id],var,status,et,eta)
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
                    ShowerShapeEvolutionPlot(can,labels,radz_hists,singlephoton_hists)
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
                    ShowerShapeEvolutionPlot(pads[-1],labels,radz_hists,singlephoton_hists)
                    plotfunc.AddHistogram(pads[-1],cuts_graphs['loose'],drawopt='l')
                    plotfunc.AddHistogram(pads[-1],cuts_graphs['tight'],drawopt='l')
                    if options.menu3 :
                        plotfunc.AddHistogram(pads[-1],cuts_graphs['menu3'],drawopt='l')
                    if options.menu4 :
                        plotfunc.AddHistogram(pads[-1],cuts_graphs['menu4'],drawopt='l')
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
                        plotfunc.MakeLegend(pads[-1],0.7,0.88,0.9,0.99,totalentries=3,option='l')
                        pads[-1].GetPrimitive('legend').SetMargin(0.4)
                        plotfunc.DrawText(pads[-1],text_lines,0.33,0.88,0.5,0.99,totalentries=3)
                    else :
                        plotfunc.DrawText(pads[-1],['','',text_lines[-1]],0.01,0.88,0.5,0.99,totalentries=1)
                        if pads[-1].GetPrimitive('legend') :
                            pads[-1].GetPrimitive('legend').Delete()

                    if eta1 == 1 :
                        plotfunc.MakeLegend(pads[-1],0.01,0.93,0.5,0.99,totalentries=2,option='l')

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
    p.add_option('--radztreename',type = 'string', default = 'output', dest = 'radztreename' ,help = 'Radiative-Z treename' )

    p.add_option('--singlephotonsignal'  ,type = 'string', default = '', dest = 'singlephotonsignal' ,help = 'Single photon Signal file' )
    p.add_option('--singlephotontreename',type = 'string', default = 'SinglePhoton', dest = 'singlephotontreename' ,help = 'Single photon treename' )

    p.add_option('--outdir',type='string',default='.',dest='outdir',help='output directory')

    options,args = p.parse_args()
    
    main(options,args)
