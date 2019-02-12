#!/usr/bin/env python

import os,sys
import ROOT
from array import array

ROOT.gROOT.SetBatch(True)

# Add base path
the_path = ('/').join(os.path.abspath(__file__).split('/')[:-2]) 

python_path = '%s/python'%(the_path)
print 'Adding %s to PYTHONPATH.'%(python_path)
sys.path.append(python_path)

genericUtils_path = '%s/genericUtils'%(the_path)
print 'Adding %s to PYTHONPATH.'%(genericUtils_path)
sys.path.append(genericUtils_path)

import python.PlotFunctions as plotfunc
import python.TAxisFunctions as taxisfunc
import python.PyAnalysisPlotting as anaplot

#-----------------------------------------------
def main(options,args) :

    plotfunc.SetupStyle()
    ROOT.gROOT.GetStyle('mystyle').SetEndErrorSize(2)

    c_u = 'c' if options.converted else 'u'
    c_uc = 'c' if options.converted else 'uc'
    un = '' if options.converted else 'un'

    can = plotfunc.RatioCanvas('Efficiency_%s_eta%s'%(c_u,options.eta),'blah',500,500,0.5)

    labels = {'mm':'Matrix Method','ee':'Electron Extrapolation','rad':'Radiative Z','comb':'Combination'}

    files = dict()
    files['mm' ] = ROOT.TFile('ID_efficiencies_matrix_method.root','READ')
    files['rad'] = ROOT.TFile('radiativeZ_data15_16_17_idtight_dep_wploose.root','READ')
    files['ee' ] = ROOT.TFile('electronExtrapolation_2015_2016_2017_closCorrFix.root','READ')

    files['sf'] = ROOT.TFile('combination_Fall2018_%s.root'%(c_u),'READ')
    eff_plots = dict()
    eff_plots['mm' ] = files['mm' ].Get('err_tot_%s_eta%s_data'%(c_u,options.eta))
    eff_plots['rad'] = files['rad'].Get('err_tot_%s_eta%s_data'%(c_u,options.eta))
    eff_plots['ee' ] = files['ee' ].Get('g_dataEff_totUnc_etabin%s_%s'%(options.eta,c_uc))

    sf_plots = dict()    
    sf_can = files['sf'].Get('%sconverted_Combined_SF'%(un))
    sf_pad = sf_can.GetPrimitive('%sconverted_Combined_SF_%d'%(un,int(options.eta)+1))
    sf_plots['mm' ]  = sf_pad.GetListOfPrimitives()[1]
    sf_plots['rad']  = sf_pad.GetListOfPrimitives()[3]
    sf_plots['ee' ]  = sf_pad.GetListOfPrimitives()[2]
    sf_plots['comb'] = sf_pad.GetListOfPrimitives()[4]

    def NewXaxis(plot) :
        add = 0
        if plot.GetName() in ['mm','ee'] :
            add = 3
        for b in range(plot.GetN()) :
            plot.GetX()[b] = b + 0.5 + add
            if issubclass(type(plot),ROOT.TGraphErrors) :
                plot.GetEX()[b] = 0.5
            if issubclass(type(plot),ROOT.TGraphAsymmErrors) :
                plot.GetEXlow()[b] = 0.5
                plot.GetEXhigh()[b] = 0.5
        return

    for k in eff_plots.keys() :
        eff_plots[k].SetName(k)
        #NewXaxis(eff_plots[k])
        eff_plots[k].SetTitle(labels[k])
        eff_plots[k].SetMarkerSize(0.8)
        eff_plots[k].SetMarkerStyle(20)
        eff_plots[k].SetLineWidth(1)
        eff_plots[k].GetHistogram().SetOption('p')

    for k in sf_plots.keys() :
        sf_plots[k].SetName(k)
        #NewXaxis(sf_plots[k])
        sf_plots[k].SetTitle(labels[k])
        sf_plots[k].SetMarkerSize(0.8)
        sf_plots[k].SetMarkerStyle(20)
        sf_plots[k].SetLineWidth(1)

    sf_plots['comb'].SetMarkerStyle(25)

    plotfunc.AddHistogram(can.GetPrimitive('pad_bot'),sf_plots['comb'],'E2')
    for method in ['mm','rad','ee'] :
        plotfunc.AddRatioManual(can,eff_plots[method],sf_plots[method],'p','p')
    sf_comb_points = ROOT.TGraph(sf_plots['comb'].GetN(),sf_plots['comb'].GetX(),sf_plots['comb'].GetY())
    sf_comb_points.SetMarkerSize(0.8)
    sf_comb_points.SetMarkerStyle(24)
    sf_comb_points.SetMarkerColor(ROOT.kBlack)
    sf_plots['comb'].SetMarkerStyle(24)
    sf_plots['comb'].SetMarkerColor(ROOT.kBlack)
    sf_plots['comb'].SetLineWidth(0)
    plotfunc.AddHistogram(can.GetPrimitive('pad_bot'),sf_comb_points,'p')

    can.GetPrimitive('pad_top').SetLogx(True)
    can.GetPrimitive('pad_bot').SetLogx(True)
    taxisfunc.SetXaxisRanges(can,9.99,1300)

    can.GetPrimitive('pad_top').SetBottomMargin(0.01)
    can.GetPrimitive('pad_bot').SetTopMargin(0.01)

    #taxisfunc.SetXaxisRanges(can,0,16)
    for p in can.GetPrimitive('pad_bot').GetListOfPrimitives() :
        if not hasattr(p,'GetXaxis') :
            continue
        p.GetXaxis().SetMoreLogLabels(True)
        p.GetXaxis().SetNoExponent(True)
        
    plotfunc.SetAxisLabels(can,'E_{T} [GeV]','data efficiency','data / MC')
    plotfunc.FormatCanvasAxes(can)

    text_lines = [plotfunc.GetAtlasInternalText('Internal')]
    text_lines += [plotfunc.GetSqrtsText(13)+', 80.0 fb^{#minus1}']
    text_lines += ['converted #gamma, '] if options.converted else ['unconverted #gamma, ']
    text_lines[-1] += ['0.0 < |#eta| < 0.6','0.6 < |#eta| < 1.37','1.52 < |#eta| < 1.81','1.81 < |#eta| < 2.37'][int(options.eta)].replace(' ','^{ }')
    plotfunc.DrawText(can,text_lines,0.45,0.05,0.80,0.28,totalentries=3)

    plotfunc.SetColors(can,[ROOT.kRed+1,ROOT.kAzure-2,ROOT.kGreen+1])
    plotfunc.MakeLegend(can,0.55,0.37,0.75,0.60,totalentries=3)
    #plotfunc.AutoFixAxes(can,ignorelegend=True)
    taxisfunc.SetYaxisRanges(can.GetPrimitive('pad_top'),0.45,1.05)
    taxisfunc.SetYaxisRanges(can.GetPrimitive('pad_bot'),0.85,1.15)

    can.GetPrimitive('pad_top').cd()
    tmp = []
    for i in can.GetPrimitive('pad_top').GetListOfPrimitives() :
        if issubclass(type(i),ROOT.TGraph) :
            tmp.append(i)
    for i in tmp :
        i.Draw('pX')

    leg = ROOT.TLegend(0.55,0.8,0.75,0.9)
    leg.SetFillStyle(0)
    leg.AddEntry(sf_plots['comb'],'Combination','pf')
    can.GetPrimitive('pad_bot').cd()
    leg.Draw()
    can.GetPrimitive('pad_bot').Update()

    options.save = True
    anaplot.doSaving(options,[can])
#     raw_input('pause')

    return

#-----------------------------------------------
if __name__ == '__main__':
    from optparse import OptionParser
    p = OptionParser()

    p.add_option('--eta'  ,type = 'string', default = '0', dest = 'eta' ,help = 'Eta bin' )
    p.add_option('--converted',action='store_true',default=False,dest='converted',help='Converted plot')
    options,args = p.parse_args()
    
    main(options,args)

