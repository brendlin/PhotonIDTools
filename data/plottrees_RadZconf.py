
treename = 'output'

bkgs = 'mc16a.Sherpa_CT10_mumugamma.root'

cuts = ['ph.convFlag==0']
# cuts = ['ph.convFlag==1']

variables = {
    'ph.pt',
    'ph.eta2',
    'ph.rhad*(ph.eta2 >= 0.8 && ph.eta2 < 1.37) + ph.rhad1*(1-(ph.eta2 >= 0.8 && ph.eta2 < 1.37))',
    'ph.e277',
    'ph.reta',
    'ph.rphi',
    'ph.weta2',
    'ph.f1',
    'ph.fside',
    'ph.wstot',
    'ph.w1',
    'ph.deltae',
    'ph.eratio',
    }

weight = 'mc_weight.pu*mc_weight.gen'

histformat = {
#     'ph.rhad*(ph.eta2 >= 0.8 && ph.eta2 < 1.37) + ph.rhad1*(1-(ph.eta2 >= 0.8 && ph.eta2 < 1.37))':['Rhad']
    }
