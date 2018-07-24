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

ROOT.gROOT.SetMacroPath('%s:%s/share'%(ROOT.gROOT.GetMacroPath(),the_path))
ROOT.gROOT.LoadMacro('CutsMacro.h')

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
        if ''.join(et_bins) :
            et_bins = list(float(a.rstrip().lstrip()) for a in et_bins)
            # Add first bin, 10 GeV
            et_bins = [10000] + et_bins
            # Remove all bins above 100 GeV
            for et in range(len(et_bins)-1,-1,-1) :
                if et_bins[et] > 100000.1 :
                    et_bins.pop(et)
            # Add 100 GeV bin threshold
            if int(et_bins[-1]) != 100000 :
                et_bins = et_bins + [100000]
        else :
            et_bins = [10000,15000,20000,25000,30000,35000,40000,45000,50000,60000,80000,100000]
        print et_bins

        args = len(et_bins)-1,array('d',list(a/1000. for a in et_bins)),len(eta_bins)-1,array('d',eta_bins)
        numerator_rz = ROOT.TH2F('numerator_%s_%s'%(status,'radz'),'Z#rightarrow^{}ll#gamma MC (fudged)',*args)
        denominator_rz = ROOT.TH2F('denominator_%s_%s'%(status,'radz'),'Z#rightarrow^{}ll#gamma MC (fudged)',*args)

        numerator_sp = ROOT.TH2F('numerator_%s_%s'%(status,'incl'),'Inclusive #gamma MC (fudged)',*args)
        denominator_sp = ROOT.TH2F('denominator_%s_%s'%(status,'incl'),'Inclusive #gamma MC (fudged)',*args)

        # Test code - for fixing up the plot cosmetics:
        if options.cosmetics :

            # TRandom3 used for plot cosmetics test code
            rand = ROOT.TRandom3()

            for et in range(len(et_bins)-1) :
                for eta in range(len(eta_bins)-1) :

                    if options.radzsignal :
                        passing = rand.Poisson(1000)
                        failing = rand.Poisson(100)
                        denominator_rz.SetBinContent(et+1,eta+1,passing+failing)
                        denominator_rz.SetBinError  (et+1,eta+1,0.01)
                        numerator_rz.SetBinContent  (et+1,eta+1,passing)
                        numerator_rz.SetBinError    (et+1,eta+1,0.01)

                    if options.singlephotonsignal :
                        passing = rand.Poisson(1000)
                        failing = rand.Poisson(100)
                        denominator_sp.SetBinContent(et+1,eta+1,passing+failing)
                        denominator_sp.SetBinError  (et+1,eta+1,0.01)
                        numerator_sp.SetBinContent  (et+1,eta+1,passing)
                        numerator_sp.SetBinError    (et+1,eta+1,0.01)

            ## End eta block
        ## End Et block

        # general setup of tight ID menu
        n_et_tight  = 1+len(confs['tight'].GetValue('CutBinEnergy_photons%s'%(status),'').split(';'))
        n_eta_tight = len(confs['tight'].GetValue('CutBinEta_photons%s'%(status),'').split(';'))
        tight_id = ROOT.photonID(n_et_tight,n_eta_tight)

        for var in variables.keys() :
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
                                                  numerator_rz)

        # New Single Photon efficiency method:
        if options.singlephotonsignal :
            print 'Evaluating Single-photon signal ID...'
            ROOT.EvaluatePhotonID_InclusivePhoton(trees_sp[keys_sp[0]],
                                                  tight_id,
                                                  status == 'Converted',
                                                  denominator_sp,
                                                  numerator_sp)

        can2d = ROOT.TCanvas('can2d_%s'%(status),'blah',600,500);
        numerator_rz.Divide(numerator_rz,denominator_rz,1,1,'B')
        numerator_sp.Divide(numerator_sp,denominator_sp,1,1,'B')
        #numerator_rz.Draw('colz')

        can = ROOT.TCanvas('can_pt_%s'%(status),'blah',600,900);
        pads = []
        hists_rz = []
        hists_sp = []

        bm = 0.08 # bottom margin
        padheight = (1-bm-0.04)/float(len(eta_bins)-2)
        for eta in range(len(eta_bins)-1) :
            if eta == 4 : continue
            eta1 = eta - 1*(eta>4)
            pads.append(ROOT.TPad('pad_eta%d','blah',0,bm*(eta1>0) + eta1*padheight,1, bm + (eta1+1)*padheight))

            current_hists = []

            # Radiative-Z: Project onto x-axis (Et)
            if options.radzsignal :
                name = numerator_rz.GetName()+'_eta%d'%(eta+1)
                numerator_rz.ProjectionX(name,eta+1,eta+1)
                hists_rz.append(ROOT.gDirectory.Get(name))
                current_hists.append(hists_rz[-1])

            # Single-photon: Project onto x-axis (Et)
            if options.singlephotonsignal :
                name = numerator_sp.GetName()+'_eta%d'%(eta+1)
                numerator_sp.ProjectionX(name,eta+1,eta+1)
                hists_sp.append(ROOT.gDirectory.Get(name))
                hists_sp[-1].SetMarkerColor(ROOT.kRed+1)
                hists_sp[-1].SetLineColor(ROOT.kRed+1)
                current_hists.append(hists_sp[-1])


            drawopt = 'pE'
            if eta1 > 0 :
                pads[-1].SetBottomMargin(0)
                pads[-1].SetTopMargin(0)
            if eta == 0 :
                pads[-1].SetTopMargin(0)
                pads[-1].SetBottomMargin(bm/float(pads[-1].GetHNDC()))                
            if eta1 == 6 :
                drawopt = 'pEX+'

            for h in current_hists :
                plotfunc.AddHistogram(pads[-1],h,drawopt)

            plotfunc.FormatCanvasAxes(pads[-1])

            for h in list(pads[-1].GetListOfPrimitives()) :
                if not issubclass(type(h),ROOT.TH1) :
                    continue

                h.SetMarkerSize(0.6)
                h.GetYaxis().SetRangeUser(.4001,0.999)
                h.GetYaxis().SetNdivisions(5,5,0)
                h.GetYaxis().SetTitleSize(16)
                h.GetYaxis().SetTitleOffset(2.3)
                h.GetXaxis().SetTitleOffset(5.2)

                if eta1 != 0 and eta1 != 6 :
                    h.GetXaxis().SetTickLength(0)

            can.cd()
            pads[-1].Draw()

            tmpy = '|#eta|^{ }#in^{ }' if eta == 0 else ''
            yaxis_label = '%s[%.2f,%.2f]'%(tmpy,eta_bins[eta],eta_bins[eta+1])
            xaxislabel  = 'p_{T}^{#gamma} [GeV]' if eta == 0 else ''
            plotfunc.SetAxisLabels(pads[-1],xaxislabel,yaxis_label)

            if eta == 1 :
                # some plot text
                text_lines = [plotfunc.GetAtlasInternalText()]
                text_lines += [plotfunc.GetSqrtsText(13)+', 13 fb^{#minus1}']
                plotfunc.DrawText(pads[-1],text_lines,0.65,0.05,0.9,0.55,totalentries=2)
            if eta == 0 :
                # Legend
                plotfunc.MakeLegend(pads[-1],.57,.47,.7,.75,option='pE',totalentries=1)

        status1 = 'Converted' if status == 'Converted' else 'Unconverted'
        text = ROOT.TLatex()
        text.SetTextFont(43)
        text.SetTextSize(18)
        text.SetTextAlign(31)
        text.DrawLatex(0.95,0.97,status1)
        can.Update()
        if not os.path.exists('%s'%(options.outdir)) :
            os.makedirs('%s'%(options.outdir))
        can.Print('%s/Efficiency_VersusPt_%s.pdf'%(options.outdir,status))

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

    p.add_option('--outdir',type='string',default='.',dest='outdir',help='output directory')

    p.add_option('--cosmetics',action='store_true',default=False,dest='cosmetics',help='For fixing plot cosmetics')

    options,args = p.parse_args()
    
    main(options,args)
