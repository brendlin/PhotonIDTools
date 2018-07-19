
def GetCutValueFromConf(conf,var,status,etbin,etabin) :
    # "conf" is a TEnv built from a typical Photon ID conf file.
    # Warning: this function assumes that the et and eta binning that you use as arguments
    # correspond to the et and eta bins in the corresponding conf.

    conf_item = '%s_photons%s'%(var,status)
    cut_value = conf.GetValue(conf_item,'')
    if not cut_value :
        return None

    # Get rid of comments wrapped with # on either side
    cut_value = ''.join(list(a if not i%2 else '' for i,a in enumerate(cut_value.split('#'))))
    try :
        cut_value = cut_value.split(';')[etbin*9 + etabin] # Et-dependent
    except IndexError :
        cut_value = cut_value.split(';')[etabin] # Et-independent
    cut_value = float(cut_value.rstrip().lstrip())

    return cut_value
