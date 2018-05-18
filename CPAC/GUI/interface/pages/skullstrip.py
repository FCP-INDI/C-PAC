#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May  7 13:47:59 2018

@author: nanditharajamani
"""

import wx
import wx.html
from ..utils.generic_class import GenericClass
from ..utils.constants import control,dtype
from ..utils.validator import CharValidator
import os
import pkg_resources as p

class SkullStripPreprocessing(wx.html.HtmlWindow):
    
    def __init__(self,parent,counter = 0):
        from urllib2 import urlopen
        wx.html.HtmlWindow.__init__(self,parent,style = wx.html.HW_SCROLLBAR_AUTO)
        self.SetStandardFonts()
        
        self.counter = counter
        self.LoadPage(p.resource_filename('CPAC', 'GUI/resources/html/anat.html'))
    
    def get_counter(self):
        return self.counter

class SkullStripOptions(wx.html.HtmlWindow):

    def __init__(self,parent,counter = 0):
        wx.ScrolledWindow.__init__(self,parent)

        self.counter = counter
        self.page = GenericClass(self, "Skull-Strip options")

        self.page.add(label ="Which function do you want to skull-strip with?",
                      control=control.CHOICE_BOX,
                      name='skullstrip_option',
                      type=dtype.LSTR,
                      comment = "Choice of using AFNI or FSL-BET to perform SkullStripping",
                      values = ["AFNI","BET"],
                      wkf_switch = True)
        self.page.set_sizer()
        parent.get_page_list().append(self)


    def get_counter(self):
        return self.counter


class AFNI_options(wx.ScrolledWindow):
    
    def __init__(self,parent,counter = 0):
        wx.ScrolledWindow.__init__(self,parent)

        
        self.counter = counter
        self.page = GenericClass(self, "AFNI options")
        
        self.page.add(label="Shrink factor",
                      control=control.TEXT_BOX,
                      name = 'shrink_factor',
                      type=dtype.LNUM,
                      comment="Set the threshold value controling the brain vs non-brain voxels\
                              default is 0.6",
                      validator=CharValidator("no-alpha"),
                      values="0.6")
                      
        self.page.add(label="Vary Shrink Factor?",
                      control=control.CHOICE_BOX,
                      name='var_shrink_fac',
                      type=dtype.LSTR,
                      comment="Vary the shrink factor at every iteration of the algorithm? this prevents the likehood of surface from getting stuck in large pools of CSF before reaching the outer surface of the brain. This is the default",
                      values=["On","Off"])
                      
        self.page.add(label="Shrink Factor Bottom Limit",
                      control=control.TEXT_BOX,
                      name='shrink_factor_bottom_lim',
                      type=dtype.LNUM,
                      comment="The shrink factor bottom limit sets the lower threshold when varying the shrink factor. Default is 0.65",
                      validator = CharValidator("no-alpha"),
                      values="0.65")
                      
        self.page.add(label="Avoid ventricles",
                      control=control.CHOICE_BOX,
                      name='avoid_vent',
                      type=dtype.LSTR,
                      comment="Avoids ventricles while skullstripping,Use this option twice for more aggressive stripping",
                      values=["On","Off"])
                      
        self.page.add(label="n-iterations",
                      control = control.TEXT_BOX,
                      name = 'n_iterations',
                      type=dtype.LNUM,
                      comment="Set the number of iterations, default is 250 and the number of iterations will depend upon the density of your mesh",
                      validator=CharValidator("no-alpha"),
                      values="250")
                      
        self.page.add(label="Pushout",
                      control=control.CHOICE_BOX,
                      name = 'pushout',
                      type = dtype.LSTR,
                      comment="While expanding, consider the voxels above and not only the voxels below",
                      values=["On","Off"])
                      
        self.page.add(label="Touchup",
                      control=control.CHOICE_BOX,
                      name = 'touchup',
                      type=dtype.LSTR,
                      comment="Perform touchup operations at the end to include areas not covered by surface expansion",
                      values=["On","Off"])
                      
        self.page.add(label = "Fill_hole",
                      control=control.TEXT_BOX,
                      name = 'fill_hole',
                      type = dtype.LNUM,
                      comment="Give the maximum number of pixels on either side of the hole that can be filled",
                      validator = CharValidator("no-alpha"),
                      values = "0")
                      
        self.page.add(label="NN_smooth",
                      control=control.TEXT_BOX,
                      name = 'NN_smooth',
                      type = dtype.LNUM,
                      comment = "Perform nearest neighbor coordinate interpolation every few iterations.Default is 72",
                      validator = CharValidator("no-alpha"),
                      values = "72")
                      
        self.page.add(label="Smooth_final",
                      control=control.TEXT_BOX,
                      name = 'smooth_final',
                      type = dtype.LNUM,
                      comment = "Perform final surface smoothing after all iterations. Default is 20",
                      validator=CharValidator("no-alpha"),
                      values = "20")
                      
        self.page.add(label="Avoid_eyes",
                      control = control.CHOICE_BOX,
                      name = 'avoid_eyes',
                      type = dtype.LSTR,
                      comment = "Avoid eyes while skull stripping,defualt is True",
                      values = ["On","Off"])
                      
        self.page.add(label="Use_edge",
                      control = control.CHOICE_BOX,
                      name = 'use_edge',
                      type = dtype.LSTR,
                      comment = "Use edge detection to reduce leakage into meninges and eyes, default is True",
                      values = ["On","Off"])
                      
        self.page.add(label = "Push_to_edge",
                      control = control.CHOICE_BOX,
                      name = 'push_to_edge',
                      type = dtype.LSTR,
                      comment = "Perform aggressive push to edge, this might cause leakage",
                      values = ["On","Off"])
                      
        self.page.add(label = "Perc_init",
                      control = control.TEXT_BOX,
                      name = 'perc_init',
                      type = dtype.LNUM,
                      comment = "Percentage of segments allowed to intersect surface. It is typically a number between 0 and 0.1, but can include negative values (which implies no testing for intersection",
                      validator=CharValidator("no-alpha"),
                      values = "0")
                      
        self.page.add(label = "Max_inter_init",
                      control = control.TEXT_BOX,
                      name = 'max_inter_init',
                      type = dtype.LNUM,
                      comment = "blur dset after spatial normalization. Recommended values are between 2 and 4",
                      validator=CharValidator("no-alpha"),
                      values = "2")
                      
        self.page.add(label = "Fac",
                      control = control.TEXT_BOX,
                      name = 'fac',
                      type = dtype.LNUM,
                      comment = "Multiply input dataset by FAC if range of values is too small",
                      validator=CharValidator("no-alpha"),
                      values = "1")

        self.page.set_sizer()
        parent.get_page_list().append(self)
            

    def get_counter(self):
        return self.counter

class BET_options(wx.ScrolledWindow):
    def __init__(self,parent,counter = 0):
        wx.ScrolledWindow.__init__(self,parent)
        
        self.counter = counter
        self.page = GenericClass(self,"BET_options")
        
        self.page.add(label="center",
                      control=control.TEXT_BOX,
                      name='center',
                      comment="identify center of gravity for the image",
                      type=dtype.LNUM,
                      validator=CharValidator("no-alpha"),
                      values="[0,0,0]")
        
        self.page.add(label="Mask_boolean",
                      control=control.CHOICE_BOX,
                      name='mask_boolean',
                      comment="mask created along with skull stripping",
                      type=dtype.LSTR,
                      values=["On","Off"])
                      
        self.page.add(label="Mesh_boolean",
                      control=control.CHOICE_BOX,
                      name='mesh_boolean',
                      comment="mesh created along with skull stripping",
                      type=dtype.LSTR,
                      values=["On","Off"])
                      
                      
        self.page.add(label="Outline",
                      control=control.CHOICE_BOX,
                      name='outline',
                      comment="create a surface outline image",
                      type=dtype.LSTR,
                      values=["On","Off"])
                      
                      
        self.page.add(label="Padding",
                      control=control.CHOICE_BOX,
                      name='padding',
                      comment="add padding to the end of the image, improving BET",
                      type=dtype.LSTR,
                      values=["On","Off"])
                      
                      
        self.page.add(label="Radius",
                      control=control.TEXT_BOX,
                      name='radius',
                      comment="integer value of head radius",
                      type=dtype.LNUM,
                      validator=CharValidator("no-alpha"),
                      values="0")
                      
                      
        self.page.add(label="Reduce_bias",
                      control=control.CHOICE_BOX,
                      name='reduce_bias',
                      comment="reduce bias and cleanup neck",
                      type=dtype.LSTR,
                      values=["On","Off"])
                      
                      
        self.page.add(label="Remove_eyes",
                      control=control.CHOICE_BOX,
                      name='remove_eyes',
                      comment="eyes and optic nerve cleanup",
                      type=dtype.LSTR,
                      values=["On","Off"])
                      
                      
        self.page.add(label="Robust",
                      control=control.CHOICE_BOX,
                      name='robust',
                      comment="robust brain center estimation",
                      type=dtype.LSTR,
                      values=["On","Off"])
                      
                      
        self.page.add(label="Skull",
                      control=control.CHOICE_BOX,
                      name='skull',
                      comment="create a skull image",
                      type=dtype.LSTR,
                      values=["On","Off"])
                      
                      
        self.page.add(label="Surfaces",
                      control=control.CHOICE_BOX,
                      name='surfaces',
                      type=dtype.LSTR,
                      comment="gets additional skull and scalp surfaces by running bet2 and betsurf",
                      values=["On","Off"])
                      
                      
        self.page.add(label="Threshold",
                      control=control.CHOICE_BOX,
                      name='threshold',
                      type=dtype.LSTR,
                      comment="apply thresholding to segmented brain image and mask",
                      values=["On","Off"])
                      
        self.page.add(label="Vertical_gradient",
                      control=control.TEXT_BOX,
                      name='vertical_gradient',
                      comment="vertical gradient in fractional intensity threshold (-1,1)",
                      type=dtype.LNUM,
                      validator=CharValidator("no-alpha"),
                      values="0.000")
        self.page.set_sizer()
        parent.get_page_list().append(self)
            

    def get_counter(self):
        return self.counter
        
                

                
