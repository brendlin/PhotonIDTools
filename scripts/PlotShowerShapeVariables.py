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

variables = {
#     'mc_weight.pu' ,
#     'mc_weight.gen',
#     'ph.convFlag'  ,
#     'ph.pt'        ,
#     'ph.eta2'      ,
#     'ph.rhad'      ,
#     'ph.e277'      ,
#     'ph.f1'        ,
    'CutHadLeakage':['ph.rhad*(0.8 < fabs(ph.eta2) && fabs(ph.eta2) < 1.37) + ph.rhad1*(!(0.8 < fabs(ph.eta2) && fabs(ph.eta2) < 1.37))',-0.03,0.07,'R_{Had}'],
    'Reta37'       :['ph.reta'      ,   0.80,    1.01,'R_{#eta}'    ],
    'Rphi33'       :['ph.rphi'      ,   0.45,     1.0,'R_{#phi}'    ],
    'weta2'        :['ph.weta2'     , 0.0065, 0.01599,'W_{#eta^{}2}'],
    'fracm'        :['ph.fside'     ,  0.001,   0.799,'f_{side}'    ],
    'wtot'         :['ph.wstot'     ,      1,   4.999,'W_{stot}'    ],
    'w1'           :['ph.w1'        ,   0.45,    0.85,'W_{s3}'      ],
    'deltae'       :['ph.deltae'    , 0.0001,   799.9,'#Delta^{}E'  ],
    'DEmaxs1'      :['ph.eratio'    , 0.6001,   1.049,'E_{ratio}'   ],
    }

#-----------------------------------------------
def main(options,args) :

    plotfunc.SetupStyle()

    files_s,trees_s,keys_s = anaplot.GetTreesFromFiles(options.signal,treename=options.treename)

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

            treevar = variables[var][0]
            bins = '(%s,%s,%s)'%(100,variables[var][1],variables[var][2])

            # Reset tobject_collector
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

                sig_hists = []

                # Tight cuts graph
                cuts_tight_gr = ROOT.TGraph(2*len(et_bins)-2,
                                         array('d',[0]*(2*len(et_bins)-2)),
                                         array('d',[0]*(2*len(et_bins)-2)))
                cuts_tight_gr.SetName('tight_cuts_%s_%s_%d'%(var,status,eta))
                cuts_tight_gr.SetTitle('Tight cuts')
                cuts_tight_gr.SetLineWidth(2)
                cuts_tight_gr.SetLineColor(ROOT.kBlue+1)

                # Loose cuts graph
                cuts_loose_gr = ROOT.TGraph(2*len(et_bins)-2,
                                         array('d',[0]*(2*len(et_bins)-2)),
                                         array('d',[0]*(2*len(et_bins)-2)))
                cuts_tight_gr.SetName('loose_cuts_%s_%s_%d'%(var,status,eta))
                cuts_loose_gr.SetTitle('Loose cuts')
                cuts_loose_gr.SetLineWidth(2)
                cuts_loose_gr.SetLineColor(ROOT.kOrange+5)

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

                    # Get the histogram
                    name = '%s_%s_%d_%d'%(var,status,et,eta)
                    arg1,arg2,arg3 = '%s>>%s%s'%(treevar,name,bins),weight,'egoff'
                    ROOT.gEnv.SetValue('Hist.Binning.1D.x','100')
                    print 'tree.Draw(\'%s\',\'%s\',\'%s\')'%(arg1,arg2,arg3)
                    tmp = ROOT.gErrorIgnoreLevel
                    ROOT.gErrorIgnoreLevel = ROOT.kFatal
                    trees_s[keys_s[0]].Draw(arg1,arg2,arg3)
                    ROOT.gErrorIgnoreLevel = tmp

                    sig_hists.append(ROOT.gDirectory.Get(name))
                    sig_hists[-1].SetTitle('signal')
                    sig_hists[-1].SetLineColor(ROOT.kBlack)
                    sig_hists[-1].SetLineWidth(2)

                    # Get the tight cut values
                    conf_item = '%s_photons%s'%(var,status)
                    cut_value = confs['tight'].GetValue(conf_item,'')
                    cut_value = ''.join(list(a if not i%2 else '' for i,a in enumerate(cut_value.split('#'))))
                    cut_value = cut_value.split(';')[et*9 + eta]
                    cut_value = float(cut_value.rstrip().lstrip())
                    cuts_tight_gr.SetPoint(et*2  ,cut_value,et  )
                    cuts_tight_gr.SetPoint(et*2+1,cut_value,et+1)

                    # now loose
                    cut_value = confs['loose'].GetValue(conf_item,'')
                    if cut_value :
                        cut_value = ''.join(list(a if not i%2 else '' for i,a in enumerate(cut_value.split('#'))))
                        cut_value = cut_value.split(';')[eta]
                        cut_value = float(cut_value.rstrip().lstrip())
                        cuts_loose_gr.SetPoint(et*2  ,cut_value,et  )
                        cuts_loose_gr.SetPoint(et*2+1,cut_value,et+1)

                # some plot text
                text_lines = [plotfunc.GetAtlasInternalText()]
                text_lines += [plotfunc.GetSqrtsText(13)+', 13 fb^{#minus1}']
                text_lines += ['%2.2f < |#eta| < %2.2f'%(eta_bins[eta],eta_bins[eta+1])]

                # Individual canvas
                if False :
                    can = ROOT.TCanvas('can_%s_%s_%s'%(var,status,eta),'blah',520,700)
                    ShowerShapeEvolutionPlot(can,labels,sig_hists)
                    plotfunc.AddHistogram(can,cuts_loose_gr,drawopt='l')
                    plotfunc.AddHistogram(can,cuts_tight_gr,drawopt='l')
                    plotfunc.SetAxisLabels(can,variables[var][3],'')
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
                    ShowerShapeEvolutionPlot(pads[-1],labels,sig_hists)
                    plotfunc.AddHistogram(pads[-1],cuts_loose_gr,drawopt='l')
                    plotfunc.AddHistogram(pads[-1],cuts_tight_gr,drawopt='l')
                    if not eta1 :
                        plotfunc.SetLeftMargin(pads[-1],0.33)
                    if eta1 :
                        plotfunc.SetLeftMargin(pads[-1],0.01)
                    if eta1 < 3 :
                        plotfunc.SetRightMargin(pads[-1],0.01)
                    can_barrel.cd() if (eta < 5) else can_endcap.cd()
                    pads[-1].Draw()
                    plotfunc.SetAxisLabels(pads[-1],variables[var][3],'')
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

        ## end status block
    return
    
#-----------------------------------------------
if __name__ == '__main__':
    from optparse import OptionParser
    p = OptionParser()
    p.add_option('--tight',type = 'string', default = '', dest = 'tight',help = 'Tight Menu')
    p.add_option('--loose' ,type = 'string', default = '', dest = 'loose' ,help = 'Loose Menu' )

    p.add_option('--signal' ,type = 'string', default = '', dest = 'signal' ,help = 'Signal file' )
    p.add_option('--treename',type = 'string', default = 'output', dest = 'treename' ,help = 'treename' )

    options,args = p.parse_args()
    
    main(options,args)
