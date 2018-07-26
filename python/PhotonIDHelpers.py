
def GetCutValuesFromConf(conf,var,status) :

    # "conf" is a TEnv built from a typical Photon ID conf file.
    # Warning: this function assumes that the et and eta binning that you use as arguments
    # correspond to the et and eta bins in the corresponding conf.

    conf_item = '%s_photons%s'%(var,status)
    cut_values = conf.GetValue(conf_item,'')
    if not cut_values :
        return []

    # Get rid of comments wrapped with # on either side
    cut_values = ''.join(list(a if not i%2 else '' for i,a in enumerate(cut_values.split('#'))))
    cut_values = cut_values.split(';')
    cut_values = list(float(a.rstrip().lstrip()) for a in cut_values)
    return cut_values


def GetCutValueFromConf(conf,var,status,etbin,etabin) :
    cut_values = GetCutValuesFromConf(conf,var,status)
    if not cut_values :
        return None
    try :
        cut_value = cut_values[etbin*9 + etabin] # Et-dependent
    except IndexError :
        cut_value = cut_values[etabin] # Et-independent

    return cut_value
