# This conf file is optimized to be used on EGAM3 Z->eegamma sherpa / data samples,
# using 43.5938 fb^-1 of data with cuts for the Z->ee+FSR(gamma) channel

treename = 'output'
fb = 43.5938

data = ''
bkgs = ''

cuts = ['ph.tight_id_dep==1','ph.isotight==1','ph.pt*0.001>15','ll.m*0.001<80','85<llg.m*0.001','llg.m*0.001<95']

variables = {
    'll.m*0.001',
    'llg.m*0.001',
    'ph.pt*0.001'
    }

weight = 'mc_weight.pu*mc_weight.gen*mc_weight.xs'
def weightscale(myfile):
    return 0.001 #need to scale this as the int.lumi. must actually be given in pb^-1

histformat = {
    'll.m*0.001':[50,40,85,'m_{ll} [GeV]'],
    'llg.m*0.001':[50,80,100,'m_{ll#gamma} [GeV]'],
    'ph.pt*0.001':[100,0,200,'p_{T}^{#gamma} [GeV]']
    }

labels = {'%Pt35_70%':'MC pt 35-70',
    '%Pt10_35%':'MC pt 10-35',
    '%Pt140%':'MC pt > 140',
    '%Pt70_140%':'MC pt 70-140'
    }

rebin = {
    'ph.pt*0.001':[0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,50,54,58,62,68,74,80,88,98,108,128,200]
}
