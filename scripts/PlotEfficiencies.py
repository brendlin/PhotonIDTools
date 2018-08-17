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

#-----------------------------------------------
def main(options,args) :

    plotfunc.SetupStyle()
    ROOT.gStyle.SetPadTickX(0)

    used_confs = []
    confs = dict()
    for confstr in ['menu1','menu2','menu3'] :
        if not getattr(options,confstr) :
            continue
        print 'Using %s for %s'%(getattr(options,confstr),confstr)
        used_confs.append(confstr)
        confs[confstr] = ROOT.TEnv(getattr(options,confstr))

    samples = []
    colors = dict()
    marker_styles = dict()
    trees = dict()
    titles = dict()

    colors['menu1'] = ROOT.kBlack
    colors['menu2'] = ROOT.kRed+1
    colors['menu3'] = ROOT.kAzure-2

    if options.radzsignal :
        samples.append('radz')
        files_rz,trees_rz,keys_rz = anaplot.GetTreesFromFiles(options.radzsignal,treename=options.radztreename)
        trees['radz'] = trees_rz[keys_rz[0]] # yeah sorry about this syntax.
        colors['radz'] = ROOT.kBlack
        titles['radz'] = 'Z#rightarrow^{}ll#gamma MC (fudged)'
        marker_styles['radz'] = 33

    if options.singlephotonsignal :
        samples.append('incl')
        files_sp,trees_sp,keys_sp = anaplot.GetTreesFromFiles(options.singlephotonsignal,treename=options.singlephotontreename)
        trees['incl'] = trees_sp[keys_sp[0]]
        colors['incl'] = ROOT.kRed+1
        titles['incl'] = 'Inclusive #gamma MC (fudged)'
        marker_styles['incl'] = 20

    if options.jetfilteredbackground :
        samples.append('jf')
        files_jf,trees_jf,keys_jf = anaplot.GetTreesFromFiles(options.jetfilteredbackground,treename=options.singlephotontreename)
        trees['jf'] = trees_jf[keys_jf[0]]
        colors['jf'] = ROOT.kAzure-2
        titles['jf'] = 'Filtered jet MC (fudged)'
        marker_styles['jf'] = 34

    for status in ['Converted','NonConverted'] :

        # Eta bins for Plots!
        eta_bins = [0,0.6,0.8,1.15,1.37,1.52,1.81,2.01,2.37]
        print eta_bins

        # Assume for now that Strips and non-strips binning is the same
        et_bins = [10000,15000,20000,25000,30000,35000,40000,45000,50000,60000,80000,100000]
        print et_bins

        args = len(et_bins)-1,array('f',list(a/1000. for a in et_bins)),len(eta_bins)-1,array('f',eta_bins)

        numerators = dict()
        denominators = dict()

        for conf in used_confs :
            numerators[conf] = dict()
            denominators[conf] = dict()
            for sample in samples :
                numerators[conf][sample]   = ROOT.TH2F('numerator_%s_%s_%s'  %(status,sample,conf),titles.get(sample),*args)
                denominators[conf][sample] = ROOT.TH2F('denominator_%s_%s_%s'%(status,sample,conf),titles.get(sample),*args)

        # Test code - for fixing up the plot cosmetics:
        if options.cosmetics :

            # TRandom3 used for plot cosmetics test code
            rand = ROOT.TRandom3()

            for et in range(len(et_bins)-1) :
                for eta in range(len(eta_bins)-1) :

                    for conf in used_confs :
                        for sample in samples :
                            passing = rand.Poisson(1000)
                            failing = rand.Poisson(100)
                            denominators[conf][sample].SetBinContent(et+1,eta+1,passing+failing)
                            denominators[conf][sample].SetBinError  (et+1,eta+1,0.01)
                            numerators[conf][sample].SetBinContent  (et+1,eta+1,passing)
                            numerators[conf][sample].SetBinError    (et+1,eta+1,0.01)

                ## End eta block
            ## End Et block

        photonIDs = dict()
        for c in used_confs :
            # general setup of tight ID menu
            # There is a general quirk with the conf files on bin thresholds vs bins:
            # Eta: nbins = nthresholds (0 is implied as a threshold)
            # Et: nbins = nthresholds + 1 (0 and inf are implied as thresholds)
            etbins_conf  = idhelpers.GetCutValuesFromConf(confs[c],'CutBinEnergy',status)
            etabins_conf =  idhelpers.GetCutValuesFromConf(confs[c],'CutBinEta',status)
            n_et  = 1 + len(etbins_conf) # see note above
            n_eta = len(etabins_conf)
            photonIDs[c] = ROOT.photonID(n_et,n_eta)
            photonIDs[c].Set_EtaBinThresholds(array('f',etabins_conf))
            if etbins_conf :
                photonIDs[c].Set_EtBinThresholds(array('f',etbins_conf))

            # Add each cut
            variables = ['CutHadLeakage','Reta37','Rphi33','weta2','fracm','wtot','w1','deltae','DEmaxs1']
            for var in variables :
                cut_values = idhelpers.GetCutValuesFromConf(confs[c],var,status)
                if not cut_values :
                    continue
                cut_values = array('f',cut_values)
                getattr(photonIDs[c],'Set_%s'%(var))(cut_values)


        # New efficiency method:
        for sample in samples :
            print 'Evaluating %s signal ID...'%(sample)
            doTruthMatchPhoton = True
            doTruthMatchFake = False
            if sample == 'jf' :
                doTruthMatchPhoton = False
                doTruthMatchFake = True

            for conf in used_confs :
                ROOT.EvaluatePhotonID(trees[sample],
                                      photonIDs[conf],
                                      status == 'Converted',
                                      denominators[conf][sample],
                                      numerators[conf][sample],
                                      None,
                                      doTruthMatchPhoton, doTruthMatchFake
                                      )

                # Divide (using Binomial errors)
                numerators[conf][sample].Divide(numerators[conf][sample],denominators[conf][sample],1,1,'B')

        ##
        ## "Compressed" plots (lots of plots on the same canvas)
        ##
        can = ROOT.TCanvas('can_pt_%s'%(status),'blah',600,900);
        pads = []
        hists_ptplot = dict()
        for conf in used_confs :
            hists_ptplot[conf] = dict()
            for sample in samples :
                hists_ptplot[conf][sample] = []

        bm = 0.08 # bottom margin
        padheight = (1-bm-0.04)/float(len(eta_bins)-2)
        miny,maxy=9e12,-9e12
        for eta in range(len(eta_bins)-1) :
            if eta == 4 : continue
            eta1 = eta - 1*(eta>4)
            pads.append(ROOT.TPad('pad_eta%d','blah',0,bm*(eta1>0) + eta1*padheight,1, bm + (eta1+1)*padheight))

            current_hists = []

            for conf in used_confs :
                for sample in samples :
                    name = numerators[conf][sample].GetName()+'_eta%d'%(eta+1)
                    numerators[conf][sample].ProjectionX(name,eta+1,eta+1)
                    hists_ptplot[conf][sample].append(ROOT.gDirectory.Get(name))

                    if len(samples) == 1 and len(used_confs) > 1 :
                        hists_ptplot[conf][sample][-1].SetTitle(conf)

                    if len(used_confs) == 1 :
                        hists_ptplot[conf][sample][-1].SetMarkerColor(colors[sample])
                        hists_ptplot[conf][sample][-1].SetLineColor(colors[sample])
                    else :
                        hists_ptplot[conf][sample][-1].SetMarkerColor(colors[conf])
                        hists_ptplot[conf][sample][-1].SetLineColor(colors[conf])
                    current_hists.append(hists_ptplot[conf][sample][-1])

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
            tmp_miny,tmp_maxy = taxisfunc.AutoFixYaxis(pads[-1])
            miny,maxy = min(miny,tmp_miny),max(maxy,tmp_maxy)

            for h in list(pads[-1].GetListOfPrimitives()) :
                if not issubclass(type(h),ROOT.TH1) :
                    continue

                h.SetMarkerSize(0.6)
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

        for pad in pads :
            miny = max(miny,0)
            taxisfunc.SetYaxisRanges(pad,miny+0.0001,maxy-0.0001)

        status1 = 'Converted' if status == 'Converted' else 'Unconverted'
        text = ROOT.TLatex()
        text.SetTextFont(43)
        text.SetTextSize(18)
        text.SetTextAlign(31)
        text.DrawLatex(0.95,0.97,status1)
        can.Update()
        if not os.path.exists('%s'%(options.outdir)) :
            os.makedirs('%s'%(options.outdir))

        outputname = 'Efficiency'
        if options.jetfilteredbackground :
            outputname = 'BkgEfficiency'
        can.Print('%s/%s_VersusPt_%s.pdf'%(options.outdir,outputname,status))
        ## End compressed plots

    ## End loop over Converted, NonConverted

    return
    
#-----------------------------------------------
if __name__ == '__main__':
    from optparse import OptionParser
    p = OptionParser()
    p.add_option('--menu1',type = 'string', default = '', dest = 'menu1',help = 'A Menu')
    p.add_option('--menu2',type = 'string', default = '', dest = 'menu2',help = 'Another Menu' )
    p.add_option('--menu3',type = 'string', default = '', dest = 'menu3',help = 'Yet Another Menu' )
    p.add_option('--names',type = 'string', default = '', dest = 'names',help = 'Menu names, comma-separated (menu1,menu2,menu3,menu4)')

    p.add_option('--radzsignal'  ,type = 'string', default = '', dest = 'radzsignal' ,help = 'Radiative-Z Signal file' )
    p.add_option('--radztreename',type = 'string', default = 'output', dest = 'radztreename' ,help = 'Radiative-Z treename' )

    p.add_option('--singlephotonsignal'  ,type = 'string', default = '', dest = 'singlephotonsignal' ,help = 'Single photon Signal file' )
    p.add_option('--singlephotontreename',type = 'string', default = 'SinglePhoton', dest = 'singlephotontreename' ,help = 'Single photon treename' )
    p.add_option('--jetfilteredbackground',type = 'string', default = '', dest = 'jetfilteredbackground' ,help = 'Jet-filtered background file' )

    p.add_option('--FixedCutLoose',action='store_true',default=False,dest='FixedCutLoose',help='apply FixedCutLoose preselection')
    p.add_option('--outdir',type='string',default='.',dest='outdir',help='output directory')

    p.add_option('--cosmetics',action='store_true',default=False,dest='cosmetics',help='For fixing plot cosmetics')

    options,args = p.parse_args()
    
    main(options,args)
