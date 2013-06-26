import sys
import os
import commands
import nibabel as nib
import nipype.pipeline.engine as pe
import CPAC.interfaces.afni.preprocess as preprocess

def check_tr():

    # imageData would have to be the image data from the funcFlow workflow, funcFlow outputspec.subject
    img = nib.load(funcFlow.outputs.outputspec.subject)
    
    # get header from image data, then extract TR information, TR is fourth item in list returned by get_zooms()
    imageHeader = img.get_header()
    imageZooms = imageHeader.get_zooms()
    header_tr = imageZooms[3]
                
    # If the TR information from header_tr (funcFlow) and convert_tr node (TR from config file)
    # do not match, prepare to update the TR information from either convert_tr or header_tr using
    # afni 3drefit, then append to func_to_mni
    if header_tr != convert_tr.outputs.tr:
        
        # Create 3drefit node
        threedrefit = pe.Node(interface=afni.preprocess(), name='update_tr')
                    
        if convert_tr.outputs.tr != None:
            workflow.connect(convert_tr, 'tr', threedrefit, 'tr')
        else if header_tr != None:
            threedrefit.inputs.tr = header_tr
        else:
            # if both are none?
            pass
        
        # Input file to be updated into 3drefit
        workflow.connect(func_to_mni, 'outputspec.mni_func', threedrefit, 'in_file')
                        
        print 'Warning: The TR information does not match between the config and subject list files.'
