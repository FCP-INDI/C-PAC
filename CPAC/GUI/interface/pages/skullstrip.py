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

class SkullStripProcessing(wx.html.HtmlWindow):

    def __init__(self,parent,counter = 0):
        wx.html.HtmlWindow.__init__(
            self,
            parent,
            style=wx.html.HW_SCROLLBAR_AUTO
        )
        self.SetStandardFonts()

        self.counter = counter
        self.LoadPage(
            p.resource_filename('CPAC', 'GUI/resources/html/anat.html')
        )

    def get_counter(self):
        return self.counter

class SkullStripOptions(wx.html.HtmlWindow):

    def __init__(self,parent,counter = 0):
        wx.ScrolledWindow.__init__(self,parent)
        import os

        self.counter = counter
        self.page = GenericClass(self, "Skull-Strip options")

        self.page.add(label="Inputs Already Skull-stripped? ",
                      control=control.CHOICE_BOX,
                      name='already_skullstripped',
                      type=dtype.LSTR,
                      comment="Disables skull-stripping on the anatomical inputs if they are already skull-stripped outside of C-PAC. Set this to On if your input images are already skull-stripped.",
                      values=["Off", "On"])

        self.page.add(label="Which tool for skull-stripping?",
                      control=control.CHOICE_BOX,
                      name='skullstrip_option',
                      type=dtype.STR,
                      comment="Choice of using AFNI or FSL-BET to perform "
                              "SkullStripping",
                      values=["AFNI","BET","AFNI & BET"],
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
                      name='skullstrip_shrink_factor',
                      type=dtype.NUM,
                      comment="Set the threshold value controlling the "
                              "brain vs non-brain voxels. Default is 0.6.",
                      validator=CharValidator("no-alpha"),
                      values="0.6")

        self.page.add(label="Vary Shrink Factor?",
                      control=control.CHOICE_BOX,
                      name='skullstrip_var_shrink_fac',
                      type=dtype.STR,
                      comment="Vary the shrink factor at every iteration of "
                              "the algorithm. This prevents the likelihood "
                              "of surface getting stuck in large pools of "
                              "CSF before reaching the outer surface of "
                              "the brain. Default is On.",
                      values=["On","Off"])

        self.page.add(label="Shrink Factor Bottom Limit",
                      control=control.TEXT_BOX,
                      name='skullstrip_shrink_factor_bot_lim',
                      type=dtype.NUM,
                      comment="The shrink factor bottom limit sets the "
                              "lower threshold when varying the shrink "
                              "factor. Default is 0.4, for when edge "
                              "detection is used (which is On by default), "
                              "otherwise the default value is 0.65.",
                      validator = CharValidator("no-alpha"),
                      values="0.4")

        self.page.add(label="Avoid ventricles",
                      control=control.CHOICE_BOX,
                      name='skullstrip_avoid_vent',
                      type=dtype.STR,
                      comment="Avoids ventricles while skullstripping.",
                      values=["On","Off"])

        self.page.add(label="Number of iterations",
                      control=control.TEXT_BOX,
                      name='skullstrip_n_iterations',
                      type=dtype.NUM,
                      comment="Set the number of iterations. Default is 250."
                              "The number of iterations should depend upon "
                              "the density of your mesh.",
                      validator=CharValidator("no-alpha"),
                      values="250")

        self.page.add(label="Pushout",
                      control=control.CHOICE_BOX,
                      name='skullstrip_pushout',
                      type=dtype.STR,
                      comment="While expanding, consider the voxels above "
                              "and not only the voxels below",
                      values=["On","Off"])

        self.page.add(label="Touchup",
                      control=control.CHOICE_BOX,
                      name='skullstrip_touchup',
                      type=dtype.STR,
                      comment="Perform touchup operations at the end to "
                              "include areas not covered by surface "
                              "expansion.",
                      values=["On","Off"])

        self.page.add(label="Fill hole option",
                      control=control.TEXT_BOX,
                      name='skullstrip_fill_hole',
                      type=dtype.NUM,
                      comment="Give the maximum number of pixels on either "
                              "side of the hole that can be filled. The "
                              "default is 10 only if 'Touchup' is On - "
                              "otherwise, the default is 0.",
                      validator = CharValidator("no-alpha"),
                      values = "10")

        self.page.add(label="NN smooth",
                      control=control.TEXT_BOX,
                      name='skullstrip_NN_smooth',
                      type=dtype.NUM,
                      comment="Perform nearest neighbor coordinate "
                              "interpolation every few iterations. Default "
                              "is 72.",
                      validator=CharValidator("no-alpha"),
                      values="72")

        self.page.add(label="Smooth final",
                      control=control.TEXT_BOX,
                      name='skullstrip_smooth_final',
                      type=dtype.NUM,
                      comment="Perform final surface smoothing after all "
                              "iterations. Default is 20.",
                      validator=CharValidator("no-alpha"),
                      values="20")

        self.page.add(label="Avoid eyes",
                      control=control.CHOICE_BOX,
                      name='skullstrip_avoid_eyes',
                      type=dtype.STR,
                      comment="Avoid eyes while skull stripping. Default is "
                              "On.",
                      values = ["On","Off"])

        self.page.add(label="Use edge",
                      control=control.CHOICE_BOX,
                      name='skullstrip_use_edge',
                      type=dtype.STR,
                      comment="Use edge detection to reduce leakage into "
                              "meninges and eyes. Default is On.",
                      values=["On","Off"])

        self.page.add(label="Fractional expansion",
                      control=control.TEXT_BOX,
                      name='skullstrip_exp_frac',
                      type=dtype.NUM,
                      comment="Speed of expansion.",
                      validator=CharValidator("no-alpha"),
                      values="0.1")

        self.page.add(label="Push to edge",
                      control=control.CHOICE_BOX,
                      name='skullstrip_push_to_edge',
                      type=dtype.STR,
                      comment="Perform aggressive push to edge. This might "
                              "cause leakage. Default is Off.",
                      values=["Off","On"])

        self.page.add(label="Use skull",
                      control=control.CHOICE_BOX,
                      name='skullstrip_use_skull',
                      type=dtype.STR,
                      comment="Use outer skull to limit expansion of surface "
                              "into the skull in case of very strong shading "
                              "artifacts. Use this only if you have leakage "
                              "into the skull.",
                      values=["Off","On"])

        self.page.add(label="Perc_int",
                      control=control.TEXT_BOX,
                      name='skullstrip_perc_int',
                      type=dtype.NUM,
                      comment="Percentage of segments allowed to intersect "
                              "surface. It is typically a number between 0 "
                              "and 0.1, but can include negative values "
                              "(which implies no testing for intersection).",
                      validator=CharValidator("no-alpha"),
                      values="0")

        self.page.add(label="Max_inter_iter",
                      control=control.TEXT_BOX,
                      name='skullstrip_max_inter_iter',
                      type=dtype.NUM,
                      comment="Number of iterations to remove intersection "
                              "problems. With each iteration, the program "
                              "automatically increases the amount of "
                              "smoothing to get rid of intersections. "
                              "Default is 4.",
                      validator=CharValidator("no-alpha"),
                      values="4")

        self.page.add(label="Fac",
                      control=control.TEXT_BOX,
                      name='skullstrip_fac',
                      type=dtype.NUM,
                      comment="Multiply input dataset by FAC if range of "
                              "values is too small.",
                      validator=CharValidator("no-alpha"),
                      values="1")

        self.page.add(label="blur_fwhm",
                      control=control.TEXT_BOX,
                      name='skullstrip_blur_fwhm',
                      type=dtype.NUM,
                      comment="Blur dataset after spatial normalization. "
                              "Recommended when you have lots of CSF in "
                              "brain and when you have protruding gyri "
                              "(finger like). If so, recommended value "
                              "range is 2-4. Otherwise, leave at 0.",
                      validator=CharValidator("no-alpha"),
                      values="0")

        self.page.set_sizer()
        parent.get_page_list().append(self)


    def get_counter(self):
        return self.counter

class BET_options(wx.ScrolledWindow):
    def __init__(self,parent,counter = 0):
        wx.ScrolledWindow.__init__(self,parent)

        self.counter = counter
        self.page = GenericClass(self,"BET_options")

        self.page.add(label="frac",
                      control=control.TEXT_BOX,
                      name = 'bet_frac',
                      type=dtype.NUM,
                      comment="Set the threshold value controling the brain vs non-brain voxels\
                              default is 0.5",
                      validator=CharValidator("no-alpha"),
                      values="0.5")

        self.page.add(label="mask",
                      control=control.CHOICE_BOX,
                      name='bet_mask_boolean',
                      comment="Mask created along with skull stripping",
                      type=dtype.STR,
                      values=["Off","On"],
                      wkf_switch = True)

        self.page.add(label="mesh",
                      control=control.CHOICE_BOX,
                      name='bet_mesh_boolean',
                      comment="Mesh created along with skull stripping",
                      type=dtype.STR,
                      values=["Off","On"],
                      wkf_switch = True)

        self.page.add(label="outline",
                      control=control.CHOICE_BOX,
                      name='bet_outline',
                      comment="Create a surface outline image",
                      type=dtype.STR,
                      values=["Off","On"],
                      wkf_switch = True)


        self.page.add(label="padding",
                      control=control.CHOICE_BOX,
                      name='bet_padding',
                      comment="Add padding to the end of the image, improving BET.Mutually exclusive with functional,reduce_bias,robust,padding,remove_eyes,surfaces",
                      type=dtype.STR,
                      values=["Off","On"],
                      wkf_switch = True)


        self.page.add(label="radius",
                      control=control.TEXT_BOX,
                      name='bet_radius',
                      comment="Integer value of head radius",
                      type=dtype.NUM,
                      validator=CharValidator("no-alpha"),
                      values="0")


        self.page.add(label="reduce_bias",
                      control=control.CHOICE_BOX,
                      name='bet_reduce_bias',
                      comment="Reduce bias and cleanup neck. Mutually exclusive with functional,reduce_bias,robust,padding,remove_eyes,surfaces",
                      type=dtype.STR,
                      values=["Off","On"],
                      wkf_switch = True)

        self.page.add(label="remove_eyes",
                      control=control.CHOICE_BOX,
                      name='bet_remove_eyes',
                      comment="Eyes and optic nerve cleanup. Mutually exclusive with functional,reduce_bias,robust,padding,remove_eyes,surfaces",
                      type=dtype.STR,
                      values=["Off","On"],
                      wkf_switch = True)


        self.page.add(label="robust",
                      control=control.CHOICE_BOX,
                      name='bet_robust',
                      comment="Robust brain center estimation. Mutually exclusive with functional,reduce_bias,robust,padding,remove_eyes,surfaces",
                      type=dtype.STR,
                      values=["Off","On"],
                      wkf_switch = True)


        self.page.add(label="skull",
                      control=control.CHOICE_BOX,
                      name='bet_skull',
                      comment="Create a skull image",
                      type=dtype.STR,
                      values=["Off","On"],
                      wkf_switch = True)


        self.page.add(label="surfaces",
                      control=control.CHOICE_BOX,
                      name='bet_surfaces',
                      type=dtype.STR,
                      comment="Gets additional skull and scalp surfaces by running bet2 and betsurf. This is mutually exclusive with reduce_bias, robust, padding, remove_eyes",
                      values=["Off","On"],
                      wkf_switch = True)


        self.page.add(label="threshold",
                      control=control.CHOICE_BOX,
                      name='bet_threshold',
                      type=dtype.STR,
                      comment="Apply thresholding to segmented brain image and mask",
                      values=["Off","On"],
                      wkf_switch = True)

        self.page.add(label="Vertical_gradient",
                      control=control.TEXT_BOX,
                      name='bet_vertical_gradient',
                      comment="Vertical gradient in fractional intensity threshold (-1,1)",
                      type=dtype.NUM,
                      validator=CharValidator("no-alpha"),
                      values="0.000")
        self.page.set_sizer()
        parent.get_page_list().append(self)


    def get_counter(self):
        return self.counter
