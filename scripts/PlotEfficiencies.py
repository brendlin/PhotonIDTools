#!/usr/bin/env python

import os,sys
import ROOT
from array import array
import math

# Add base path
the_path = ('/').join(os.path.abspath(__file__).split('/')[:-2]) 

python_path = '%s/python'%(the_path)
print 'Adding %s to PYTHONPATH.'%(python_path)
sys.path.append(python_path)

genericUtils_path = '%s/genericUtils'%(the_path)
print 'Adding %s to PYTHONPATH.'%(genericUtils_path)
sys.path.append(genericUtils_path)

ROOT.gROOT.SetBatch(True)

import python.PlotFunctions as plotfunc
import python.TAxisFunctions as taxisfunc
import python.PyAnalysisPlotting as anaplot
from python.ShowerShapeEvolutionPlotter import ShowerShapeEvolutionPlot
import PhotonIDHelpers as idhelpers

variables = {
#     'ph.e277'      ,
#     'ph.f1'        ,
    'CutHadLeakage':['ph.rhad*(0.8 < fabs(ph.eta2) && fabs(ph.eta2) < 1.37) + ph.rhad1*(!(0.8 < fabs(ph.eta2) && fabs(ph.eta2) < 1.37))',
                     'y_Rhad*(0.8 < fabs(y_eta_cl_s2) && fabs(y_eta_cl_s2) < 1.37) + y_Rhad1*(!(0.8 < fabs(y_eta_cl_s2) && fabs(y_eta_cl_s2) < 1.37))',
                     '<'],
    'Reta37'       :['ph.reta'      ,'y_Reta'  ,'>'],
    'Rphi33'       :['ph.rphi'      ,'y_Rphi'  ,'>'],
    'weta2'        :['ph.weta2'     ,'y_weta2' ,'<'],
    'fracm'        :['ph.fside'     ,'y_fracs1','<'],
    'wtot'         :['ph.wstot'     ,'y_wtots1','<'],
    'w1'           :['ph.w1'        ,'y_weta1' ,'<'],
    'deltae'       :['ph.deltae'    ,'y_deltae','<'],
    'DEmaxs1'      :['ph.eratio'    ,'y_Eratio','>'],
    }

#-----------------------------------------------
def main(options,args) :

    plotfunc.SetupStyle()
    ROOT.gStyle.SetPadTickX(0)

    confs = dict()
    for confstr in ['tight'] :
        print 'Using %s for %s'%(getattr(options,confstr),confstr)
        confs[confstr] = ROOT.TEnv(getattr(options,confstr))

    files_rz,trees_rz,keys_rz = anaplot.GetTreesFromFiles(options.radzsignal,treename=options.radztreename)
    files_sp,trees_sp,keys_sp = anaplot.GetTreesFromFiles(options.singlephotonsignal,treename=options.singlephotontreename)

    for status in ['Converted','NonConverted'] :

        # Assume for now that Strips and non-strips binning is the same
        eta_bins = confs['tight'].GetValue('CutBinEta_photons%s'%(status),'').split(';')
        eta_bins = list(float(a.rstrip().lstrip()) for a in eta_bins)
        eta_bins = [0] + eta_bins[:-1]
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
        numerator_rz = ROOT.TH2F('numerator_%s_%s'%(status,'radz'),'Z#rightarrow^{}ll#gamma MC (fudged)',*args)
        numerator_sp = ROOT.TH2F('numerator_%s_%s'%(status,'incl'),'Inclusive photons MC (fudged)',*args)
        denominator_rz = ROOT.TH2F('denominator_%s_%s'%(status,'radz'),'Z#rightarrow^{}ll#gamma MC (fudged)',*args)
        denominator_sp = ROOT.TH2F('denominator_%s_%s'%(status,'incl'),'Inclusive photons MC (fudged)',*args)

        for et in range(len(et_bins)-1) :
            for eta in range(len(eta_bins)-1) :
                phasespace_cuts = []
                id_cuts = []

                phasespace_cuts.append('ph.pt > %2.1f'%(et_bins[et]))
                phasespace_cuts.append('ph.pt < %2.1f'%(et_bins[et+1]))
                phasespace_cuts.append('fabs(ph.eta2) > %2.2f'%(eta_bins[eta]))
                phasespace_cuts.append('fabs(ph.eta2) < %2.2f'%(eta_bins[eta+1]))
                if status == 'Converted' :
                    phasespace_cuts.append('ph.convFlag > 0')
                else :
                    phasespace_cuts.append('ph.convFlag == 0')

                for var in variables.keys() :
                    # Get the cut values
                    cut_value = idhelpers.GetCutValueFromConf(confs['tight'],var,status,et,eta)
                    if cut_value == None :
                        continue
                    id_cuts.append('%s %s %s'%(variables[var][0],variables[var][2],cut_value))

                # options for the histogram-getter function
                options.limits = dict()
                options.limits['fabs(ph.eta2)'] = [1,0,5]

                # First, Radiative-Z
                # Get the denominator number
                weight_radz = 'mc_weight.pu*mc_weight.gen'
                weight = (weight_radz+'*(%s)'%(' && '.join(phasespace_cuts))).lstrip('*')
                hist = anaplot.GetVariableHistsFromTrees(trees_rz,keys_rz,'fabs(ph.eta2)',weight,options)
                if hist :
                    hist = hist[0]
                    Ntotal = hist.Integral(0,hist.GetNbinsX()+1)
                    Dtotal = math.sqrt(sum(list(hist.GetSumw2())))
                    denominator_rz.SetBinContent(et+1,eta+1,Ntotal)
                    denominator_rz.SetBinError  (et+1,eta+1,Dtotal)

                # Get the numerator number
                weight = (weight_radz+'*(%s)'%(' && '.join(phasespace_cuts + id_cuts))).lstrip('*')
                hist = anaplot.GetVariableHistsFromTrees(trees_rz,keys_rz,'fabs(ph.eta2)',weight,options)
                if hist :
                    hist = hist[0]
                    Npass = hist.Integral(0,hist.GetNbinsX()+1)
                    Dpass = math.sqrt(sum(list(hist.GetSumw2())))
                    numerator_rz.SetBinContent(et+1,eta+1,Npass)
                    numerator_rz.SetBinError  (et+1,eta+1,Dpass)


            ## End eta block
        ## End Et block

        can2d = ROOT.TCanvas('can2d_%s'%(status),'blah',600,500);
        numerator_rz.Divide(numerator_rz,denominator_rz,1,1,'B')
        numerator_rz.Draw('colz')

        can = ROOT.TCanvas('can_pt_%s'%(status),'blah',600,700);
        pads = []
        hists = []
        bm = 0.08 # bottom margin
        padheight = (1-bm-0.04)/float(len(eta_bins)-2)
        for eta in range(len(eta_bins)-1) :
            if eta == 4 : continue
            eta1 = eta - 1*(eta>4)
            pads.append(ROOT.TPad('pad_eta%d','blah',0,bm*(eta1>0) + eta1*padheight,1, bm + (eta1+1)*padheight))
            name = numerator_rz.GetName()+'_eta%d'%(eta+1)
            numerator_rz.ProjectionX(name,eta+1,eta+1)
            hists.append(ROOT.gDirectory.Get(name))

            drawopt = 'pE'
            if eta1 > 0 :
                pads[-1].SetBottomMargin(0)
                pads[-1].SetTopMargin(0)
            if eta == 0 :
                pads[-1].SetTopMargin(0)
                pads[-1].SetBottomMargin(bm/float(pads[-1].GetHNDC()))                
            if eta1 == 6 :
                drawopt = 'pEX+'

            pads[-1].cd()
            hists[-1].Draw(drawopt)
            hists[-1].SetDrawOption(drawopt)
            plotfunc.FormatCanvasAxes(pads[-1])

            hists[-1].SetMarkerSize(0.6)
            hists[-1].GetYaxis().SetRangeUser(.4001,0.999)
            hists[-1].GetYaxis().SetNdivisions(5,2,0)
            hists[-1].GetYaxis().SetTitleSize(16)
            hists[-1].GetYaxis().SetTitleOffset(2.3)
            hists[-1].GetXaxis().SetTitleOffset(5.2)
            can.cd()
            pads[-1].Draw()

            if eta1 != 0 and eta1 != 6 :
                hists[-1].GetXaxis().SetTickLength(0)

            tmpy = '|#eta|^{ }#in^{ }' if eta == 0 else ''
            yaxis_label = '%s[%.2f,%.2f]'%(tmpy,eta_bins[eta],eta_bins[eta+1])
            xaxislabel  = 'p_{T}^{#gamma} [GeV]' if eta == 0 else ''
            plotfunc.SetAxisLabels(pads[-1],xaxislabel,yaxis_label)

            if eta == 1 :
                # some plot text
                text_lines = [plotfunc.GetAtlasInternalText()]
                text_lines += [plotfunc.GetSqrtsText(13)+', 13 fb^{#minus1}']
                plotfunc.DrawText(pads[-1],text_lines,0.65,0.1,0.9,0.55,totalentries=2)
            if eta == 0 :
                # Legend
                plotfunc.MakeLegend(pads[-1],.65,.6,.7,.7,option='pE',totalentries=1)

        status1 = 'Converted' if status == 'Converted' else 'Unconverted'
        text_lines2 = ['%s photons'%(status)]
        plotfunc.DrawText(can,text_lines2,0.71,0.96,0.91,0.99,totalentries=1)
        can.Update()
        can.Print('Efficiency_PtDependent_%s.pdf'%(status))

    return
    
#-----------------------------------------------
if __name__ == '__main__':
    from optparse import OptionParser
    p = OptionParser()
    p.add_option('--tight',type = 'string', default = '', dest = 'tight',help = 'Tight Menu')
    p.add_option('--loose' ,type = 'string', default = '', dest = 'loose' ,help = 'Loose Menu' )
    p.add_option('--menu3' ,type = 'string', default = '', dest = 'menu3' ,help = 'Another Menu' )

    p.add_option('--radzsignal'  ,type = 'string', default = '', dest = 'radzsignal' ,help = 'Radiative-Z Signal file' )
    p.add_option('--radztreename',type = 'string', default = 'output', dest = 'radztreename' ,help = 'Radiative-Z treename' )

    p.add_option('--singlephotonsignal'  ,type = 'string', default = '', dest = 'singlephotonsignal' ,help = 'Single photon Signal file' )
    p.add_option('--singlephotontreename',type = 'string', default = 'SinglePhoton', dest = 'singlephotontreename' ,help = 'Single photon treename' )

    options,args = p.parse_args()
    
    main(options,args)

