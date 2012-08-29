def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('CPAC', parent_package, top_path)
    
    #Not sure why this isn't done in other python packages
    #this might not be stable?
    import CPAC
    for p in CPAC.__all__:
        config.add_subpackage(p)
    
    config.add_subpackage('interfaces')
    config.add_subpackage('interfaces/afni')
    
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())