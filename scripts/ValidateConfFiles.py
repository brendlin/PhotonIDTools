#!/usr/bin/env python

##
## Use this script to validate whether Loose, Medium and Tight photon
## conf files are subsets are one another, whenever validating any
## new cut-based ID.
##

import ROOT

smaller_is_tighter = {'CutHadLeakage': True,
                      'CutF3'        : True,
                      'Reta37'       : False,
                      'Rphi33'       : False,
                      'weta2'        : True,
                      'f1'           : False,
                      'deltae'       : True,
                      'DEmaxs1'      : False,
                      'wtot'         : True,
                      'fracm'        : True,
                      'w1'           : True,
                      }

class color:
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'

#-----------------------------------------------
def main(options,args) :
    
    confs = dict()
    values_not_checked = dict()
    
    print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    for confstr in ['tighter','looser'] :
        confs[confstr] = ROOT.TEnv(getattr(options,confstr))
        values_not_checked[confstr] = list(a.GetName() for a in confs[confstr].GetTable())

        print '~~ %s (%s)'%(confstr,getattr(options,confstr))
        #confs[confstr].Print()
    print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        
    # Check if tight offline is tighter than tight online
    should_be_identical = ['CutBinEta_photonsNonConverted',
                           'CutBinEtaStrips_photonsNonConverted',
                           'CutBinEnergyStrips_photonsNonConverted',
                           'CutBinEnergy_photonsNonConverted',
                           'e277_photonsNonConverted', # Determines if other cuts are applied too
                           ]

    for i in range(len(should_be_identical)) :
        should_be_identical.append(should_be_identical[i].replace('NonConverted','Converted'))

    for i in should_be_identical :
        offline = confs['tighter'].GetValue(i,'')
        online  = confs['looser' ].GetValue(i,'')

        if (not online and not offline) :
            try :
                values_not_checked['tighter'].pop(values_not_checked['tighter'].index(i))
            except ValueError :
                pass
            try :
                values_not_checked['looser' ].pop(values_not_checked['looser' ].index(i))
            except ValueError :
                pass
            continue

        if (offline and not online) :
            print 'ERROR: %s missing for online.'%(i)
            values_not_checked['tighter'].pop(values_not_checked['tighter'].index(i))
            continue

        if (online and not offline) :
            print 'ERROR: %s missing for offline.'%(i)
            values_not_checked['looser' ].pop(values_not_checked['looser' ].index(i))
            continue

        l_offline = list(float(a) for a in offline.split(';'))
        l_online  = list(float(a) for a in online .split(';'))

        if offline == online :
            values_not_checked['tighter'].pop(values_not_checked['tighter'].index(i))
            values_not_checked['looser' ].pop(values_not_checked['looser' ].index(i))
        elif l_offline == l_online :
            values_not_checked['tighter'].pop(values_not_checked['tighter'].index(i))
            values_not_checked['looser' ].pop(values_not_checked['looser' ].index(i))            
        else :
            print 'ERROR: Check %s'%(i)
            print 'tighter:',offline
            print 'looser:',online
            import sys; sys.exit()

    cuts = ['CutHadLeakage_photonsNonConverted',
            'CutF3_photonsNonConverted',
            'Reta37_photonsNonConverted',
            'Rphi33_photonsNonConverted',
            'weta2_photonsNonConverted',
            'f1_photonsNonConverted',
            'deltae_photonsNonConverted',
            'DEmaxs1_photonsNonConverted',
            'wtot_photonsNonConverted',
            'fracm_photonsNonConverted',
            'w1_photonsNonConverted',
            ]

    for i in range(len(cuts)) :
        cuts.append(cuts[i].replace('NonConverted','Converted'))

    for i in cuts :

        cut_shortname = i.replace('_photons','').replace('NonConverted','').replace('Converted','')
        conv_unconv = 'Unconv' if 'NonConverted' in i else 'Conv'
        cut_shortname_bold = color.BOLD + cut_shortname + ' (' + conv_unconv + ')' + color.END
        
        offline = confs['tighter'].GetValue(i,'')
        online  = confs['looser' ].GetValue(i,'')

        if not offline and not online :
            continue

        if not online :
            print 'INFO: %s is not applied in Online.'%(cut_shortname_bold)
            values_not_checked['tighter'].pop(values_not_checked['tighter'].index(i))            
            continue

        l_offline = list(float(a) for a in offline.split(';'))
        l_online  = list(float(a) for a in online .split(';'))

        if len(l_offline) != len(l_online) :
            print 'ERROR: One conf has the wrong length!\n',i,'Offline:',l_offline,'\n',i,'Online :',l_online
            import sys; sys.exit()

        bad_indices = []

        for j in range(len(l_offline)) :
            if smaller_is_tighter[cut_shortname] and l_offline[j] > l_online[j] :
                bad_indices.append(j)
            if (not smaller_is_tighter[cut_shortname]) and l_offline[j] < l_online[j] :
                bad_indices.append(j)

        str_offline = ', '.join((color.BOLD + ('%.6f'%(a)).rstrip('0').rjust(9) + color.END if (i in bad_indices) else ('%.6f'%(a)).rstrip('0').rjust(9) ) for i,a in enumerate(l_offline))
        str_online  = ', '.join((color.BOLD + ('%.6f'%(a)).rstrip('0').rjust(9) + color.END if (i in bad_indices) else ('%.6f'%(a)).rstrip('0').rjust(9) ) for i,a in enumerate(l_online ))
        
        str_larger_smaller = 'larger' if smaller_is_tighter[cut_shortname] else 'smaller'
        str_i = i.ljust(34)

        if len(bad_indices) :
            print 'ERROR: %s Offline is looser (%s) than online!\n%s Offline: %s\n%s Online : %s'%(cut_shortname_bold,str_larger_smaller,str_i,str_offline,str_i,str_online)
        else :
            print 'INFO: %s is ok!\n%s Offline: %s\n%s Online : %s'%(cut_shortname_bold,str_i,str_offline,str_i,str_online)

        values_not_checked['tighter'].pop(values_not_checked['tighter'].index(i))
        values_not_checked['looser' ].pop(values_not_checked['looser' ].index(i))

    if len(values_not_checked['tighter']) or len(values_not_checked['looser']) :
        print 'ERROR: Offline Not checked:', values_not_checked['tighter']
        print 'ERROR: Online  Not checked:', values_not_checked['looser' ]

    return

#-----------------------------------------------
if __name__ == '__main__':
    from optparse import OptionParser
    p = OptionParser()
    p.add_option('--tighter',type = 'string', default = '', dest = 'tighter',help = 'Menu that is supposed to be tighter')
    p.add_option('--looser' ,type = 'string', default = '', dest = 'looser' ,help = 'Menu that is supposed to be looser' )
    options,args = p.parse_args()
    
    main(options,args)
