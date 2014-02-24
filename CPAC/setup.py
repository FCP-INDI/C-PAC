def configuration(parent_package='', top_path=None):
    import os
    from numpy.distutils.misc_util import Configuration
    
    config = Configuration('CPAC', parent_package, top_path)
    
    # This adds all the first level sub-packages
    _,dirnames,_ = os.walk('CPAC').next()
    for dirname in dirnames:
        _,_,filenames = os.walk('CPAC/%s' % dirname).next()
        if '__init__.py' in filenames:
            config.add_subpackage(dirname)
    
    config.add_subpackage('interfaces/afni')
    config.add_subpackage('GUI/interface')
    config.add_subpackage('GUI/interface/pages')
    config.add_subpackage('GUI/interface/windows')
    config.add_subpackage('GUI/interface/utils')
    
    # Data
    config.add_data_dir('GUI/resources')
    config.add_data_dir('resources')
    
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
