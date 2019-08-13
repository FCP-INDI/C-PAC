

def add_afni_prefix(tpattern):
    if tpattern:
        if ".txt" in tpattern:
            tpattern = "@{0}".format(tpattern)
    return tpattern