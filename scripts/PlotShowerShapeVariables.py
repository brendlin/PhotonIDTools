#!/usr/bin/env python

import os,sys
import ROOT
from array import array

# Add base path
the_path = ('/').join(os.path.abspath(__file__).split('/')[:-2]) 
genericUtils_path = '%s/genericUtils'%(the_path)
print 'Adding %s to PYTHONPATH.'%(genericUtils_path)
sys.path.append(genericUtils_path)
ROOT.gROOT.SetMacroPath(genericUtils_path)

ROOT.gROOT.SetBatch(True)

import python.PlotFunctions as plotfunc
import python.TAxisFunctions as taxisfunc
import python.PyAnalysisPlotting as anaplot
from python.ShowerShapeEvolutionPlotter import ShowerShapeEvolutionPlot

weight_radz = 'mc_weight.pu*mc_weight.gen'
weight_singlephoton = 'mcTotWeightNoPU_PIDuse'

variables = {
#     'ph.e277'      ,
#     'ph.f1'        ,
    'CutHadLeakage':['ph.rhad*(0.8 < fabs(ph.eta2) && fabs(ph.eta2) < 1.37) + ph.rhad1*(!(0.8 < fabs(ph.eta2) && fabs(ph.eta2) < 1.37))',
                     'y_Rhad*(0.8 < fabs(y_eta_cl_s2) && fabs(y_eta_cl_s2) < 1.37) + y_Rhad1*(!(0.8 < fabs(y_eta_cl_s2) && fabs(y_eta_cl_s2) < 1.37))',-0.03,0.07,'R_{Had}'],
    'Reta37'       :['ph.reta'      ,'y_Reta'  ,   0.80,    1.01,'R_{#eta}'    ],
    'Rphi33'       :['ph.rphi'      ,'y_Rphi'  ,   0.45,     1.0,'R_{#phi}'    ],
    'weta2'        :['ph.weta2'     ,'y_weta2' , 0.0065, 0.01599,'W_{#eta^{}2}'],
    'fracm'        :['ph.fside'     ,'y_fracs1',  0.001,   0.799,'f_{side}'    ],
    'wtot'         :['ph.wstot'     ,'y_wtots1',      1,   4.999,'W_{stot}'    ],
    'w1'           :['ph.w1'        ,'y_weta1' ,   0.45,    0.85,'W_{s3}'      ],
    'deltae'       :['ph.deltae'    ,'y_deltae', 0.0001,   799.9,'#Delta^{}E'  ],
    'DEmaxs1'      :['ph.eratio'    ,'y_Eratio', 0.6001,   1.049,'E_{ratio}'   ],
    }

#-----------------------------------------------
def main(options,args) :

    plotfunc.SetupStyle()

    confs = dict()
    for confstr in ['tight','loose'] :
        print 'Using %s for %s'%(getattr(options,confstr),confstr)
        confs[confstr] = ROOT.TEnv(getattr(options,confstr))

    for status in ['Converted','NonConverted'] :

        print 'Making plots for %s'%(status)

        # Assume for now that Strips and non-strips binning is the same
        eta_bins = confs['tight'].GetValue('CutBinEta_photons%s'%(status),'').split(';')
        eta_bins = list(float(a.rstrip().lstrip()) for a in eta_bins)
        eta_bins = [0] + eta_bins
        print eta_bins

        # Assume for now that Strips and non-strips binning is the same
        et_bins = confs['tight'].GetValue('CutBinEnergy_photons%s'%(status),'').split(';')
        et_bins = list(float(a.rstrip().lstrip()) for a in et_bins)
        et_bins = [0] + et_bins
        tmp = []
        for i in et_bins :
            if i < 100000.1 :
                tmp.append(i)
        et_bins = tmp
        print et_bins

        for var in variables.keys() :

            # Reset PlotFunctions tobject_collector
            del plotfunc.tobject_collector[:]

            can_barrel = ROOT.TCanvas('can_%s_%s_barrel'%(var,status),'blah',int(1464),700)
            can_endcap = ROOT.TCanvas('can_%s_%s_endcap'%(var,status),'blah',int(1464),700)
            pads = []

            files_rz,trees_rz,keys_rz = anaplot.GetTreesFromFiles(options.radzsignal,treename=options.radztreename)
            files_sp,trees_sp,keys_sp = anaplot.GetTreesFromFiles(options.singlephotonsignal,treename=options.singlephotontreename)

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
                    tmp_gr.SetLineWidth(2)
                    return tmp_gr

                cuts_graphs = dict()

                # Tight cuts graph
                cuts_graphs['tight'] = GetCutsGraph('Tight cuts')
                cuts_graphs['tight'].SetLineColor(ROOT.kBlue+1)

                # Loose cuts graph
                cuts_graphs['loose'] = GetCutsGraph('Loose cuts')
                cuts_graphs['loose'].SetLineColor(ROOT.kOrange+5)

                for et in range(len(et_bins)-1) :

                    cuts = []
                    cuts.append('ph.pt > %2.1f'%(et_bins[et]))
                    cuts.append('ph.pt < %2.1f'%(et_bins[et+1]))
                    cuts.append('ph.eta2 > %2.2f'%(eta_bins[eta]))
                    cuts.append('ph.eta2 < %2.2f'%(eta_bins[eta+1]))
                    if status == 'Converted' :
                        cuts.append('ph.convFlag > 0')
                    else :
                        cuts.append('ph.convFlag == 0')

                    weight = (weight_radz+'*(%s)'%(' && '.join(cuts))).lstrip('*')
                    treevar = variables[var][0]

                    # Needed to pass options to histogram-getter function
                    options.limits = dict()
                    options.limits[treevar] = [100,variables[var][2],variables[var][3]]

                    # Get the histogram (Rad-Z)
                    hist = anaplot.GetVariableHistsFromTrees(trees_rz,keys_rz,treevar,weight,options)
                    if hist :
                        radz_hists.append(hist[0])
                        radz_hists[-1].SetTitle('Radiative-Z')
                        radz_hists[-1].SetLineColor(ROOT.kBlack)
                        radz_hists[-1].SetLineWidth(2)

                    # Switch to SinglePhoton inputs
                    for i in range(len(cuts)) :
                        cuts[i] = cuts[i].replace('ph.pt'      ,'y_pt*1000.' )
                        cuts[i] = cuts[i].replace('ph.eta2'    ,'y_eta_cl_s2')
                        cuts[i] = cuts[i].replace('ph.convFlag','y_convType' )
                    weight = (weight_singlephoton+'*(%s)'%(' && '.join(cuts))).lstrip('*')
                    treevar = variables[var][1]
                    options.limits[treevar] = [100,variables[var][2],variables[var][3]]

                    # Get the histogram (Single-photon)
                    hist = anaplot.GetVariableHistsFromTrees(trees_sp,keys_sp,treevar,weight,options)
                    if hist :
                        singlephoton_hists.append(hist[0])
                        singlephoton_hists[-1].SetTitle('Single photon')
                        singlephoton_hists[-1].SetLineColor(ROOT.kRed+1)
                        singlephoton_hists[-1].SetLineWidth(2)

                    # Get the cut values
                    for id in ['tight','loose'] :
                        conf_item = '%s_photons%s'%(var,status)
                        cut_value = confs[id].GetValue(conf_item,'')
                        if not cut_value :
                            continue
                        cut_value = ''.join(list(a if not i%2 else '' for i,a in enumerate(cut_value.split('#'))))
                        try :
                            cut_value = cut_value.split(';')[et*9 + eta] # Et-dependent
                        except IndexError :
                            cut_value = cut_value.split(';')[eta] # Et-independent
                        cut_value = float(cut_value.rstrip().lstrip())
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
                    plotfunc.SetAxisLabels(can,variables[var][4],'')
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
                    if not eta1 :
                        plotfunc.SetLeftMargin(pads[-1],0.33)
                    if eta1 :
                        plotfunc.SetLeftMargin(pads[-1],0.01)
                    if eta1 < 3 :
                        plotfunc.SetRightMargin(pads[-1],0.01)
                    can_barrel.cd() if (eta < 5) else can_endcap.cd()
                    pads[-1].Draw()
                    plotfunc.SetAxisLabels(pads[-1],variables[var][4],'')
                    taxisfunc.SetNdivisions(pads[-1],5,5,0)
                    for i in pads[-1].GetListOfPrimitives() :
                        if hasattr(i,'SetLineWidth') :
                            i.SetLineWidth(1)
                        if hasattr(i,'GetYaxis') :
                            i.GetYaxis().SetLabelSize(19)
                            i.GetXaxis().SetLabelSize(19)
                    if not eta1 :
                        plotfunc.MakeLegend(pads[-1],0.7,0.88,0.9,0.99,totalentries=3,option='l')
                        plotfunc.DrawText(pads[-1],text_lines,0.33,0.88,0.5,0.99,totalentries=3)
                    else :
                        plotfunc.DrawText(pads[-1],['','',text_lines[-1]],0.01,0.88,0.5,0.99,totalentries=1)
                        if pads[-1].GetPrimitive('legend') :
                            pads[-1].GetPrimitive('legend').Delete()

                ## end barrel / EC block
            ## end eta block

            if not os.path.exists(status) :
                os.makedirs(status)
            can_barrel.Print('%s/barrel_%s_%s.pdf'%(status,var,status))
            can_endcap.Print('%s/endcap_%s_%s.pdf'%(status,var,status))

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

    p.add_option('--radzsignal'  ,type = 'string', default = '', dest = 'radzsignal' ,help = 'Radiative-Z Signal file' )
    p.add_option('--radztreename',type = 'string', default = 'output', dest = 'radztreename' ,help = 'Radiative-Z treename' )

    p.add_option('--singlephotonsignal'  ,type = 'string', default = '', dest = 'singlephotonsignal' ,help = 'Single photon Signal file' )
    p.add_option('--singlephotontreename',type = 'string', default = 'SinglePhoton', dest = 'singlephotontreename' ,help = 'Single photon treename' )

    options,args = p.parse_args()
    
    main(options,args)
