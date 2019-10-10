
def add_afni_prefix(tpattern):
    if tpattern:
        if ".txt" in tpattern:
            tpattern = "@{0}".format(tpattern)
    return tpattern


def nullify(value, function=None):
    from traits.trait_base import Undefined
    if value is None:
        return Undefined
    if function:
        return function(value)
    return value