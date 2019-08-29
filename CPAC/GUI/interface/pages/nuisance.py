# -*- coding: utf-8 -*-

import copy
import wx
import wx.html
import wx.lib.newevent
import wx.lib.scrolledpanel as scrolled
from wx.lib.mixins.listctrl import CheckListCtrlMixin, ListCtrlAutoWidthMixin
from ..utils.generic_class import GenericClass, Control
from ..utils.constants import control, dtype
from ..utils.validator import CharValidator

import re
import pkg_resources as p
from collections import defaultdict

def find(lst, lmbd):
    return next((t for t in lst if lmbd(t)), None)

def findi(lst, lmbd):
    return next((i for i, t in enumerate(lst) if lmbd(t)), None)

EditorOkEvent, EVT_EDITOR_OK = wx.lib.newevent.NewEvent()
EditorCancelEvent, EVT_EDITOR_CANCEL = wx.lib.newevent.NewEvent()

EnableEvent, EVT_ENABLE = wx.lib.newevent.NewEvent()
DisableEvent, EVT_DISABLE = wx.lib.newevent.NewEvent()

def ctrl_enable(ctrl):
    try:
        ctrl.SetEditable(True)
    except:
        pass
    ctrl.Enable()
    ctrl.SetBackgroundColour(wx.WHITE)
    
    evt = EnableEvent()
    wx.PostEvent(ctrl, evt)

def ctrl_disable(ctrl):
    try:
        ctrl.SetEditable(False)
    except:
        pass
    ctrl.Disable()
    color = wx.SystemSettings.GetColour(wx.SYS_COLOUR_BACKGROUND)
    ctrl.SetBackgroundColour(color)
    
    evt = DisableEvent()
    wx.PostEvent(ctrl, evt)

def dd():
    return defaultdict(dd)

class CheckListCtrl(wx.ListCtrl, CheckListCtrlMixin, ListCtrlAutoWidthMixin):
    def __init__(self, parent, **kwargs):
        wx.ListCtrl.__init__(self, parent, wx.ID_ANY, style=wx.LC_REPORT | wx.SUNKEN_BORDER, **kwargs)
        CheckListCtrlMixin.__init__(self)
        ListCtrlAutoWidthMixin.__init__(self)


class Nuisance(wx.html.HtmlWindow):

    def __init__(self, parent, counter  = 0):
        from urllib2 import urlopen
        wx.html.HtmlWindow.__init__(self, parent, style= wx.html.HW_SCROLLBAR_AUTO)
        self.SetStandardFonts()

        self.counter = counter
        self.LoadPage(p.resource_filename('CPAC', 'GUI/resources/html/nuisance.html'))

    def get_counter(self):
        return self.counter


class FilteringSettings(wx.ScrolledWindow):
    pass


class Scrubbing(wx.ScrolledWindow):
    pass


def selectors_repr(selectors):

    representation = ''
    renaming = {
        'WhiteMatter': 'WM',
        'CerebrospinalFluid': 'CSF',
        'GreyMatter': 'GM',
        'GlobalSignal': 'GR',
        'Motion': 'Mot',
    }

    regressors = []
    for reg in ['aCompCor', 'tCompCor',
                'WhiteMatter', 'CerebrospinalFluid', 'GreyMatter',
                'GlobalSignal', 'Motion']:

        if reg not in selectors:
            continue

        renamed = reg if reg not in renaming else renaming[reg]
        regressor = selectors[reg] or {}

        terms = [renamed]
        if regressor.get('include_squared'):
            terms += ["{}{}".format(renamed, '²')]
        if regressor.get('include_delayed'):
            terms += ["{}{}".format(renamed, 'ₜ₋₁')]
        if regressor.get('include_delayed_squared'):
            terms += ["{}{}".format(renamed, 'ₜ₋₁²')]

        regressor_terms = ' + '.join(terms)

        if regressor.get('tissues'):
            regressor_terms += ' ('
            tissues = [renaming[t] for t in regressor.get('tissues')]
            regressor_terms += ' + '.join(tissues) + ')'

        if 'summary' in regressor:
            if type(regressor['summary']) == dict:
                regressor_terms += ' ' + regressor['summary']['method']
                if regressor['summary']['method'] in ['DetrendPC', 'PC']:
                    regressor_terms += ' {}'.format(regressor['summary']['components'])
            else:
                regressor_terms += ' ' + regressor['summary']

        regressors += [regressor_terms]

    representation += ", ".join(regressors) + "\n"

    renaming = {
        'SpikeRegression': 'Spike Regr'
    }

    censor = selectors.get('Censor')
    if censor:
        method = censor['method']
        method = renaming[method] if method in renaming else method

        threshs = []
        for thresh in censor['thresholds']:
            threshs += ["{0} > {1}".format(thresh['type'], thresh['value'])]

        representation += "{} ({})".format(method, ", ".join(threshs))

    bandpass = selectors.get('Bandpass')
    polynomial = selectors.get('PolyOrt')

    if bandpass or polynomial:
        representation += "\n"

    if bandpass:
        representation += "Bandpass "

        top = bandpass.get('top_frequency')
        bot = bandpass.get('bottom_frequency')

        if top and bot:
            representation += "{}Hz–{}Hz".format(top, bot)
        elif top:
            representation += "{}Hz–∞Hz".format(top)
        elif bot:
            representation += "0.0Hz–{}Hz".format(bot)

    if polynomial:
        if bandpass:
            representation += ", "
        representation += "Polynomial Reg degree {}".format(polynomial['degree'])

    return representation


selector_renaming = {
    'GlobalSignal': 'Global Signal',
    'Motion': 'Motion',
    'CerebrospinalFluid': 'Cerebrospinal Fluid',
    'WhiteMatter': 'White Matter',
    'GreyMatter': 'Grey Matter',
    'aCompCor': 'aCompCor',
    'tCompCor': 'tCompCor',
    'Custom': 'Custom',
    'PolyOrt': 'Poly Regression',
    'Bandpass': 'Bandpass',
    'Censor': 'Censor',
}

selector_renaming_inverse = \
    dict(zip(*zip(*selector_renaming.items())[::-1]))

selector_summary_method_renaming = {
    'PC': 'Principal Components',
    'DetrendPC': 'Detrended Principal Components',
    'Mean': 'Mean',
    'NormMean': 'Norm Mean',
    'DetrendNormMean': 'Detrended Norm Mean',
}

selector_summary_method_renaming_inverse = \
    dict(zip(*zip(*selector_summary_method_renaming.items())[::-1]))

selector_derivatives_renaming = {
    'include_delayed': 'Include Delayed',
    'include_squared': 'Include Squared',
    'include_delayed_squared': 'Include Delayed Squared',
}

selector_derivatives_renaming_inverse = \
    dict(zip(*zip(*selector_derivatives_renaming.items())[::-1]))

selector_censor_method_renaming = {
    'Interpolate': 'Interpolate',
    'SpikeRegression': 'Spike Regression',
    'Kill': 'Kill',
    'Zero': 'Zero',
}

selector_censor_method_renaming_inverse = \
    dict(zip(*zip(*selector_censor_method_renaming.items())[::-1]))


class NuisanceRegressionRegressorEditor(wx.Frame):

    def __init__(self, parent, selectors_index=None):
        wx.Frame.__init__(self, parent,
                          title="Edit Selectors",
                          size=(650, 400))

        self.selectors_index = selectors_index
        self.selectors = {}
        if selectors_index is not None:
            self.selectors = copy.deepcopy(parent.regressor_selectors[selectors_index])

        root = wx.Panel(self)
        sizer = wx.BoxSizer(wx.VERTICAL)
        root.SetSizer(sizer)


        panel = wx.Panel(root)
        panel_sizer = wx.BoxSizer(wx.HORIZONTAL)
        panel.SetSizer(panel_sizer)

        # Render tree
        self.tree = wx.TreeCtrl(panel,
                                style=wx.TR_HIDE_ROOT | wx.TR_FULL_ROW_HIGHLIGHT | wx.TR_SINGLE)
        self.Bind(wx.EVT_TREE_SEL_CHANGED, self.scroll_selector, self.tree)

        # Render Editor
        editor = wx.Panel(panel)
        editor_sizer = wx.BoxSizer(wx.VERTICAL)
        editor.SetSizer(editor_sizer)

        self.editor = scrolled.ScrolledPanel(editor, style=wx.VSCROLL | wx.BORDER_SUNKEN)
        self.editor.SetBackgroundColour(wx.WHITE)
        self.editor.SetupScrolling(scroll_x=False, scroll_y=True)
        self.editor.SetBackgroundColour(wx.WHITE)
        self.editor.SetSizer(wx.GridBagSizer(vgap=5, hgap=5))

        editor_sizer.Add(self.editor, 1, wx.EXPAND)

        # Render panel w tree + editor
        panel_border = wx.TOP | wx.LEFT | wx.RIGHT

        panel_sizer.Add(self.tree, 1, wx.EXPAND | panel_border, border=5)
        panel_sizer.Add(editor, 2, wx.EXPAND | panel_border, border=5)

        sizer.Add(panel, 1, wx.EXPAND)

        # Render buttons
        button_panel = wx.Panel(root)
        button_sizer = wx.BoxSizer(wx.HORIZONTAL)
        button_panel.SetSizer(button_sizer)

        button_ok = wx.Button(button_panel, -1, 'OK', size=(90, 30))
        button_ok.Bind(wx.EVT_BUTTON, self.ok_click)
        button_sizer.Add(button_ok, 0, wx.ALIGN_CENTER)

        button_cancel = wx.Button(button_panel, -1, 'Cancel', size=(90, 30))
        button_cancel.Bind(wx.EVT_BUTTON, self.cancel_click)
        button_sizer.Add(button_cancel, 0, wx.ALIGN_CENTER)

        sizer.Add(button_panel, 0, wx.ALIGN_CENTER | wx.ALL, border=5)

        self.render_tree()
        self.render_selectors()

    def render_selectors_title(self, title, main=False):
        sizer = self.editor.GetSizer()
        font = wx.Font(14 if main else 12, wx.DEFAULT, wx.NORMAL, wx.FONTWEIGHT_BOLD)
        description = wx.StaticText(self.editor, wx.ID_ANY, title)
        description.SetFont(font)
        description.SetBackgroundColour(wx.RED)
        sizer.Add(
            description,
            pos=(sizer.Rows, 0),
            span=(1, 2),
            flag=wx.EXPAND | wx.ALL, border=5
        )
        if main:
            sizer.Add(
                wx.StaticLine(self.editor, style=wx.LI_HORIZONTAL),
                pos=(sizer.Rows, 0),
                span=(1, 2),
                flag=wx.EXPAND | wx.LEFT| wx.RIGHT| wx.BOTTOM, border=5
            )
        return description

    def add_to_new_row(self, label, control,
                       label_flag=wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT | wx.ALL,
                       control_flag=wx.ALIGN_CENTER_VERTICAL | wx.EXPAND | wx.ALL):
        sizer = self.editor.GetSizer()
        rows = sizer.Rows
        sizer.Add(label, pos=(rows, 0), flag=label_flag, border=5)
        sizer.Add(control, pos=(rows, 1), flag=control_flag, border=5)
        return rows

    def add_enabled_row(self, label, is_enabled, fields_to_control):
        enabled_label = wx.StaticText(self.editor, label=label)
        enabled_control = wx.CheckBox(self.editor)
        enabled_control.SetValue(is_enabled)
        self.add_to_new_row(enabled_label, enabled_control)

        def iterate(value, d):
            if isinstance(d, dict):
                d = d.values()
            for c in d:
                if c == value:
                    continue
                if isinstance(c, dict):
                    iterate(value, c)
                    continue
                if isinstance(c, list):
                    iterate(value, c)
                    continue
                if value.GetValue():
                    ctrl_enable(c)
                else:
                    ctrl_disable(c)

        def enabled_update(event=None):
            iterate(enabled_control, fields_to_control)

        enabled_control.Bind(wx.EVT_CHECKBOX, enabled_update)
        enabled_control.UpdateState = enabled_update

        return enabled_control

    def add_summary(self, selector, controls):
        label = wx.StaticText(self.editor, label='Summary')

        summary = wx.Panel(self.editor)
        summary_sizer = wx.BoxSizer(wx.VERTICAL)
        summary.SetSizer(summary_sizer)

        if not selector:
            selector = {}
        selector_summary = selector.get('summary', {'method': 'Mean'})
        if type(selector_summary) != dict:
            selector_summary = {'method': selector_summary}
        selector_summary_method = selector_summary.get('method', 'Mean')

        control = wx.ComboBox(summary, choices=selector_summary_method_renaming.values(), style=wx.CB_READONLY)
        control.SetValue(selector_summary_method_renaming[selector_summary_method])
        summary_sizer.Add(control, flag=wx.EXPAND | wx.ALL)

        self.add_to_new_row(label, summary)
        controls['summary']['method'] = control

        def control_update(event=None):
            if selector_summary_method_renaming_inverse.get(control.GetValue()) in ['DetrendPC', 'PC']:
                components_control = wx.TextCtrl(summary, value='1')
                summary_sizer.Add(components_control, flag=wx.EXPAND | wx.ALL)
                controls['summary']['components'] = components_control
            else:
                for l, c in list(controls['summary'].items()):
                    if l == 'method':
                        continue
                    c.Destroy()
                    del controls['summary'][l]
            
            # summary.Layout()
            summary.Parent.Layout()

        control.Bind(wx.EVT_COMBOBOX, control_update)
        control_update()

        return controls

    def add_tissue_parameters(self, selector):
        controls = {}

        extraction_resolution = None
        erode = None
        if type(selector) == dict:
            extraction_resolution = selector.get('extraction_resolution')
            erode = selector.get('erode')

        extraction_label = wx.StaticText(self.editor, label='Extraction Resolution (mm)')
        extraction_control = wx.TextCtrl(self.editor, value='2')
        if extraction_resolution is not None:
            extraction_control.SetValue(str(extraction_resolution))
        self.add_to_new_row(extraction_label, extraction_control)
        controls['extraction'] = extraction_control

        erosion_label = wx.StaticText(self.editor, label='Erode Mask')
        erosion_control = wx.CheckBox(self.editor)
        if erode is not None:
            erosion_control.SetValue(erode)
        self.add_to_new_row(erosion_label, erosion_control)
        controls['erode'] = erosion_control

        return controls

    def add_derivatives(self, selector):
        if not selector:
            selector = {}

        controls = {}
        for did, derivative in selector_derivatives_renaming.items():
            label = wx.StaticText(self.editor, label=derivative)
            control = wx.CheckBox(self.editor)
            if selector.get(did):
                control.SetValue(True)
            self.add_to_new_row(label, control)
            controls[did] = control
        return controls

    def render_selectors(self):

        sizer = self.editor.GetSizer()

        self.titles = []

        self.data_controls = dd()

        self.titles += [self.render_selectors_title('Regressors', main=True)]
        sizer.AddGrowableCol(1)
        

        self.titles += [self.render_selectors_title('Motion')]

        has_motion = 'Motion' in self.selectors

        self.data_controls['Motion']['enabled'] = self.add_enabled_row('Enabled', has_motion, self.data_controls['Motion'])
        self.data_controls['Motion'].update(self.add_derivatives(self.selectors.get('Motion')))
        self.data_controls['Motion']['enabled'].UpdateState()


        self.titles += [self.render_selectors_title('Cerebrospinal Fluid')]

        has_cerebrospinalfluid = 'CerebrospinalFluid' in self.selectors

        self.data_controls['CerebrospinalFluid']['enabled'] = self.add_enabled_row('Enabled', has_cerebrospinalfluid, self.data_controls['CerebrospinalFluid'])
        self.data_controls['CerebrospinalFluid'].update(self.add_summary(self.selectors.get('CerebrospinalFluid'), self.data_controls['CerebrospinalFluid']))
        self.data_controls['CerebrospinalFluid'].update(self.add_tissue_parameters(self.selectors.get('CerebrospinalFluid')))
        self.data_controls['CerebrospinalFluid'].update(self.add_derivatives(self.selectors.get('CerebrospinalFluid')))
        self.data_controls['CerebrospinalFluid']['enabled'].UpdateState()


        self.titles += [self.render_selectors_title('White Matter')]

        has_whitematter = 'WhiteMatter' in self.selectors

        self.data_controls['WhiteMatter']['enabled'] = self.add_enabled_row('Enabled', has_whitematter, self.data_controls['WhiteMatter'])
        self.data_controls['WhiteMatter'].update(self.add_summary(self.selectors.get('WhiteMatter'), self.data_controls['WhiteMatter']))
        self.data_controls['WhiteMatter'].update(self.add_tissue_parameters(self.selectors.get('WhiteMatter')))
        self.data_controls['WhiteMatter'].update(self.add_derivatives(self.selectors.get('WhiteMatter')))
        self.data_controls['WhiteMatter']['enabled'].UpdateState()


        self.titles += [self.render_selectors_title('Grey Matter')]

        has_greymatter = 'GreyMatter' in self.selectors

        self.data_controls['GreyMatter']['enabled'] = self.add_enabled_row('Enabled', has_greymatter, self.data_controls['GreyMatter'])
        self.data_controls['GreyMatter'].update(self.add_summary(self.selectors.get('GreyMatter'), self.data_controls['GreyMatter']))
        self.data_controls['GreyMatter'].update(self.add_tissue_parameters(self.selectors.get('GreyMatter')))
        self.data_controls['GreyMatter'].update(self.add_derivatives(self.selectors.get('GreyMatter')))
        self.data_controls['GreyMatter']['enabled'].UpdateState()


        self.titles += [self.render_selectors_title('Global Signal')]

        has_globalsignal = 'GlobalSignal' in self.selectors

        self.data_controls['GlobalSignal']['enabled'] = self.add_enabled_row('Enabled', has_globalsignal, self.data_controls['GlobalSignal'])
        self.data_controls['GlobalSignal'].update(self.add_summary(self.selectors.get('GlobalSignal'), self.data_controls['GlobalSignal']))
        self.data_controls['GlobalSignal'].update(self.add_derivatives(self.selectors.get('GlobalSignal')))
        self.data_controls['GlobalSignal']['enabled'].UpdateState()


        self.titles += [self.render_selectors_title('aCompCor')]

        has_acompcor = 'aCompCor' in self.selectors

        self.data_controls['aCompCor']['enabled'] = self.add_enabled_row('Enabled', has_acompcor, self.data_controls['aCompCor'])

        acompcor_tissues_label = wx.StaticText(self.editor, label='Tissues') 
        acompcor_tissues_control = CheckListCtrl(self.editor, size=(-1, 80))
        acompcor_tissues_control.InsertColumn(0, 'Tissues')
        acompcor_tissues_control.InsertStringItem(0, 'White Matter')
        acompcor_tissues_control.InsertStringItem(1, 'Cerebrospinal Fluid')
        self.add_to_new_row(acompcor_tissues_label, acompcor_tissues_control)
        self.data_controls['aCompCor']['tissues'] = acompcor_tissues_control

        acompcor_components_label = wx.StaticText(self.editor, label='Principal Components')
        acompcor_components_control = wx.TextCtrl(self.editor, value='1')
        acompcor_components_control.SetValue(str(self.selectors.get('aCompCor', {}).get('summary', {}).get('components', 1)))
        self.add_to_new_row(acompcor_components_label, acompcor_components_control)
        self.data_controls['aCompCor']['summary']['components'] = acompcor_components_control
        
        self.data_controls['aCompCor'].update(self.add_tissue_parameters(self.selectors.get('aCompCor')))
        self.data_controls['aCompCor'].update(self.add_derivatives(self.selectors.get('aCompCor')))
        self.data_controls['aCompCor']['enabled'].UpdateState()


        self.titles += [self.render_selectors_title('tCompCor')]

        has_tcompcor = 'tCompCor' in self.selectors

        self.data_controls['tCompCor']['enabled'] = self.add_enabled_row('Enabled', has_tcompcor, self.data_controls['tCompCor'])

        tcompcor_components_label = wx.StaticText(self.editor, label='Principal Components')
        tcompcor_components_control = wx.TextCtrl(self.editor, value='1')
        tcompcor_components_control.SetValue(str(self.selectors.get('tCompCor', {}).get('summary', {}).get('components', 1)))
        self.add_to_new_row(tcompcor_components_label, tcompcor_components_control)
        self.data_controls['tCompCor']['summary']['components'] = tcompcor_components_control

        tcompcor_threshold_label = wx.StaticText(self.editor, label='Threshold')
        tcompcor_threshold_control = wx.TextCtrl(self.editor, value='1.0')
        tcompcor_threshold_control.SetValue(str(self.selectors.get('tCompCor', {}).get('threshold', 1.0)))
        self.add_to_new_row(tcompcor_threshold_label, tcompcor_threshold_control)
        self.data_controls['tCompCor']['threshold'] = tcompcor_threshold_control

        tcompcor_byslice_label = wx.StaticText(self.editor, label='Compute by Slice')
        tcompcor_byslice_control = wx.CheckBox(self.editor)
        tcompcor_byslice_control.SetValue(self.selectors.get('tCompCor', {}).get('by_slice', False))
        self.add_to_new_row(tcompcor_byslice_label, tcompcor_byslice_control)
        self.data_controls['tCompCor']['by_slice'] = tcompcor_byslice_control

        self.data_controls['tCompCor'].update(self.add_derivatives(self.selectors.get('tCompCor')))
        self.data_controls['tCompCor']['enabled'].UpdateState()


        self.titles += [self.render_selectors_title('Custom')]

        has_custom = 'Custom' in self.selectors

        self.data_controls['Custom']['enabled'] = self.add_enabled_row('Enabled', has_custom, self.data_controls['Custom'])

        custom_components_label = wx.StaticText(self.editor, label='Regressor file')
        custom_components_control = wx.TextCtrl(self.editor, value='1')
        custom_components_control.SetValue(str(self.selectors.get('Custom', {}).get('file', '')))
        self.add_to_new_row(custom_components_label, custom_components_control)
        self.data_controls['Custom']['file'] = custom_components_control

        self.data_controls['Custom']['enabled'].UpdateState()


        self.titles += [self.render_selectors_title('Detrending and Filtering', main=True)]

        self.titles += [self.render_selectors_title('Poly Regression')]

        has_poly = 'PolyOrt' in self.selectors

        self.data_controls['PolyOrt']['enabled'] = self.add_enabled_row('Enabled', has_poly, self.data_controls['PolyOrt'])

        poly_degree_label = wx.StaticText(self.editor, label='Degree')
        poly_degree_control = wx.TextCtrl(self.editor, value='0')
        poly_degree_control.SetValue(str(self.selectors.get('PolyOrt', {}).get('degree', 2)))
        self.add_to_new_row(poly_degree_label, poly_degree_control)
        self.data_controls['PolyOrt']['degree'] = poly_degree_control

        self.data_controls['PolyOrt']['enabled'].UpdateState()


        self.titles += [self.render_selectors_title('Bandpass')]

        has_bandpass = 'Bandpass' in self.selectors

        self.data_controls['Bandpass']['enabled'] = self.add_enabled_row('Enabled', has_bandpass, self.data_controls['Bandpass'])

        bandpass_freq_bottom_label = wx.StaticText(self.editor, label='Bottom Frequency')
        bandpass_freq_bottom_control = wx.TextCtrl(self.editor, value='0.0')
        bandpass_freq_bottom_control.SetValue(str(self.selectors.get('Bandpass', {}).get('bottom_frequency', 0.0)))
        self.add_to_new_row(bandpass_freq_bottom_label, bandpass_freq_bottom_control)
        self.data_controls['Bandpass']['bottom_frequency'] = bandpass_freq_bottom_control

        bandpass_freq_top_label = wx.StaticText(self.editor, label='Top Frequency')
        bandpass_freq_top_control = wx.TextCtrl(self.editor, value='9999.9')
        bandpass_freq_top_control.SetValue(str(self.selectors.get('Bandpass', {}).get('top_frequency', 9999.9)))
        self.add_to_new_row(bandpass_freq_top_label, bandpass_freq_top_control)
        self.data_controls['Bandpass']['top_frequency'] = bandpass_freq_top_control
        
        self.data_controls['Bandpass']['enabled'].UpdateState()


        self.titles += [self.render_selectors_title('Censoring', main=True)]

        censoring = self.selectors.get('Censor', {})
        censoring_method = censoring.get('method', [])
        censoring_thresholds = censoring.get('thresholds', [])

        has_censoring = 'Censor' in self.selectors

        self.data_controls['Censor']['enabled'] = self.add_enabled_row('Enabled', has_censoring, self.data_controls['Censor'])

        censoring_method_label = wx.StaticText(self.editor, label='Method')
        censoring_method_control = wx.ComboBox(self.editor, choices=selector_censor_method_renaming.values(), style=wx.CB_READONLY)
        if censoring_method:
            censoring_method_control.SetValue(selector_censor_method_renaming[censoring_method])
        self.add_to_new_row(censoring_method_label, censoring_method_control)
        self.data_controls['Censor']['method'] = censoring_method_control

        censoring_prev_tr_label = wx.StaticText(self.editor, label='Previous TRs to censor')
        censoring_prev_tr_control = wx.TextCtrl(self.editor, value='0')
        censoring_prev_tr_control.SetValue(str(self.selectors.get('Censor', {}).get('number_of_previous_trs_to_censor', 0)))
        self.add_to_new_row(censoring_prev_tr_label, censoring_prev_tr_control)
        self.data_controls['Censor']['number_of_previous_trs_to_censor'] = censoring_prev_tr_control

        censoring_subs_tr_label = wx.StaticText(self.editor, label='Subsequent TRs to censor')
        censoring_subs_tr_control = wx.TextCtrl(self.editor, value='0')
        censoring_subs_tr_control.SetValue(str(self.selectors.get('Censor', {}).get('number_of_subsequent_trs_to_censor', 0)))
        self.add_to_new_row(censoring_subs_tr_label, censoring_subs_tr_control)
        self.data_controls['Censor']['number_of_subsequent_trs_to_censor'] = censoring_subs_tr_control

        # Framewise displacement J control
        censoring_thresholds_fdj = find(censoring_thresholds, lambda t:  t.get('type') == 'FD_J') or {}
        has_censoring_thresholds_fdj = bool(censoring_thresholds_fdj)

        censoring_thresholds_fdj_label = wx.StaticText(self.editor, label='Framewise Displacement (Jenkinson)')
        censoring_thresholds_fdj_control = wx.CheckBox(self.editor)
        censoring_thresholds_fdj_control.SetValue(has_censoring_thresholds_fdj)
        self.add_to_new_row(censoring_thresholds_fdj_label, censoring_thresholds_fdj_control)
        self.data_controls['Censor']['threshold']['fdj']['enabled'] = censoring_thresholds_fdj_control

        censoring_thresholds_fdj_value_label = wx.StaticText(self.editor, label='Threshold')
        censoring_thresholds_fdj_value_control = wx.TextCtrl(self.editor)
        censoring_thresholds_fdj_value_control.SetValue(str(censoring_thresholds_fdj.get('value', 0.0)))
        self.add_to_new_row(censoring_thresholds_fdj_value_label, censoring_thresholds_fdj_value_control)
        self.data_controls['Censor']['threshold']['fdj']['value'] = censoring_thresholds_fdj_value_control

        def censoring_thresholds_fdj_update(event=None):
            if censoring_thresholds_fdj_control.GetValue():
                ctrl_enable(censoring_thresholds_fdj_value_control)
            else:
                ctrl_disable(censoring_thresholds_fdj_value_control)

        censoring_thresholds_fdj_control.Bind(wx.EVT_CHECKBOX, censoring_thresholds_fdj_update)
        censoring_thresholds_fdj_control.Bind(EVT_ENABLE, censoring_thresholds_fdj_update)
        censoring_thresholds_fdj_update()

        # Framewise displacement P control
        censoring_thresholds_fdp = find(censoring_thresholds, lambda t:  t.get('type') == 'FD_P') or {}
        has_censoring_thresholds_fdp = bool(censoring_thresholds_fdp)

        censoring_thresholds_fdp_label = wx.StaticText(self.editor, label='Framewise Displacement (Power)')
        censoring_thresholds_fdp_control = wx.CheckBox(self.editor)
        censoring_thresholds_fdp_control.SetValue(has_censoring_thresholds_fdp)
        self.add_to_new_row(censoring_thresholds_fdp_label, censoring_thresholds_fdp_control)
        self.data_controls['Censor']['threshold']['fdp']['enabled'] = censoring_thresholds_fdp_control

        censoring_thresholds_fdp_value_label = wx.StaticText(self.editor, label='Threshold')
        censoring_thresholds_fdp_value_control = wx.TextCtrl(self.editor)
        censoring_thresholds_fdp_value_control.SetValue(str(censoring_thresholds_fdp.get('value', 0.0)))
        self.add_to_new_row(censoring_thresholds_fdp_value_label, censoring_thresholds_fdp_value_control)
        self.data_controls['Censor']['threshold']['fdp']['value'] = censoring_thresholds_fdp_value_control

        def censoring_thresholds_fdp_update(event=None):
            if censoring_thresholds_fdp_control.GetValue():
                ctrl_enable(censoring_thresholds_fdp_value_control)
            else:
                ctrl_disable(censoring_thresholds_fdp_value_control)

        censoring_thresholds_fdp_control.Bind(wx.EVT_CHECKBOX, censoring_thresholds_fdp_update)
        censoring_thresholds_fdp_control.Bind(EVT_ENABLE, censoring_thresholds_fdp_update)
        censoring_thresholds_fdp_update()

        # DVARS control
        censoring_thresholds_dvars = find(censoring_thresholds, lambda t:  t.get('type') == 'DVARS') or {}
        has_censoring_thresholds_dvars = bool(censoring_thresholds_dvars)

        censoring_thresholds_dvars_label = wx.StaticText(self.editor, label='DVARS')
        censoring_thresholds_dvars_control = wx.CheckBox(self.editor)
        censoring_thresholds_dvars_control.SetValue(has_censoring_thresholds_dvars)
        self.add_to_new_row(censoring_thresholds_dvars_label, censoring_thresholds_dvars_control)
        self.data_controls['Censor']['threshold']['dvars']['enabled'] = censoring_thresholds_dvars_control

        censoring_thresholds_dvars_value_label = wx.StaticText(self.editor, label='Threshold')
        censoring_thresholds_dvars_value_control = wx.TextCtrl(self.editor)
        censoring_thresholds_dvars_value_control.SetValue(str(censoring_thresholds_dvars.get('value', 0.0)))
        self.add_to_new_row(censoring_thresholds_dvars_value_label, censoring_thresholds_dvars_value_control)
        self.data_controls['Censor']['threshold']['dvars']['value'] = censoring_thresholds_dvars_value_control

        def censoring_thresholds_dvars_update(event=None):
            if censoring_thresholds_dvars_control.GetValue():
                ctrl_enable(censoring_thresholds_dvars_value_control)
            else:
                ctrl_disable(censoring_thresholds_dvars_value_control)

        censoring_thresholds_dvars_control.Bind(wx.EVT_CHECKBOX, censoring_thresholds_dvars_update)
        censoring_thresholds_dvars_control.Bind(EVT_ENABLE, censoring_thresholds_dvars_update)
        censoring_thresholds_dvars_update()

        self.data_controls['Censor']['enabled'].UpdateState()

    def scroll_selector(self, event):
        item_name = str(self.tree.GetItemText(event.GetItem()))
        for title in self.titles:
            if item_name == title.GetLabelText():
                _, sppy = self.editor.GetScrollPixelsPerUnit()
                sy = self.editor.GetScrollPos(wx.VERTICAL) * sppy
                _, y = title.GetPosition().Get()
                self.editor.Scroll(0, int((sy + y) / sppy) - 1)

    def render_tree(self):
        self.tree.DeleteAllItems()

        root = self.tree.AddRoot('Selectors')
        root_regressors = self.tree.AppendItem(root, 'Regressors')
        root_detrending = self.tree.AppendItem(root, 'Detrending')
        root_censoring = self.tree.AppendItem(root, 'Censoring')

        valid_regressors = [
            'Motion',
            'CerebrospinalFluid',
            'WhiteMatter',
            'GreyMatter',
            'GlobalSignal',
            'aCompCor',
            'tCompCor',
            'Custom',
        ]

        for k in valid_regressors:
            self.tree.AppendItem(root_regressors, selector_renaming[k])

        valid_detrending = [
            'PolyOrt',
            'Bandpass',
        ]

        for k in valid_detrending:
            self.tree.AppendItem(root_detrending, selector_renaming[k])

        self.tree.ExpandAll()

    def compile_selector_derivatives(self, selector):
        return {
            'include_delayed': bool(selector['include_delayed'].GetValue()),
            'include_delayed_squared': bool(selector['include_delayed_squared'].GetValue()),
            'include_squared': bool(selector['include_squared'].GetValue()),
        }

    def compile_selector_tissue_parameters(self, selector):
        return {
            'extraction_resolution': float(selector['extraction'].GetValue().replace('mm', '')),
            'erode_mask': bool(selector['erode'].GetValue()),
        }

    def compile_selector_summary(self, selector):
        method = 'PC'
        if 'method' in selector['summary']:
            method = str(selector['summary']['method'].GetValue())
        params = {
            'summary': {
                'method': method
            }
        }
        if method in ['DetrendPC', 'PC']:
            params['summary']['components'] = int(selector['summary']['components'].GetValue())
        return params

    def parse_threshold(self, threshold):
        threshold_sd = \
            re.match(r"([0-9]*\.*[0-9]*)\s*SD", str(threshold))
        if not threshold_sd:
            threshold = float(threshold)
        return threshold

    def compile_selector(self):
        selector = {}

        if self.data_controls['Motion']['enabled'].GetValue():
            selector['Motion'] = {}
            selector['Motion'].update(self.compile_selector_derivatives(self.data_controls['Motion']))

        if self.data_controls['CerebrospinalFluid']['enabled'].GetValue():
            selector['CerebrospinalFluid'] = {}
            selector['CerebrospinalFluid'].update(self.compile_selector_tissue_parameters(self.data_controls['CerebrospinalFluid']))
            selector['CerebrospinalFluid'].update(self.compile_selector_summary(self.data_controls['CerebrospinalFluid']))
            selector['CerebrospinalFluid'].update(self.compile_selector_derivatives(self.data_controls['CerebrospinalFluid']))
            
        if self.data_controls['WhiteMatter']['enabled'].GetValue():
            selector['WhiteMatter'] = {}
            selector['WhiteMatter'].update(self.compile_selector_tissue_parameters(self.data_controls['WhiteMatter']))
            selector['WhiteMatter'].update(self.compile_selector_summary(self.data_controls['WhiteMatter']))
            selector['WhiteMatter'].update(self.compile_selector_derivatives(self.data_controls['WhiteMatter']))
            
        if self.data_controls['GreyMatter']['enabled'].GetValue():
            selector['GreyMatter'] = {}
            selector['GreyMatter'].update(self.compile_selector_tissue_parameters(self.data_controls['GreyMatter']))
            selector['GreyMatter'].update(self.compile_selector_summary(self.data_controls['GreyMatter']))
            selector['GreyMatter'].update(self.compile_selector_derivatives(self.data_controls['GreyMatter']))
            
        if self.data_controls['GlobalSignal']['enabled'].GetValue():
            selector['GlobalSignal'] = {}
            selector['GlobalSignal'].update(self.compile_selector_summary(self.data_controls['GlobalSignal']))
            selector['GlobalSignal'].update(self.compile_selector_derivatives(self.data_controls['GlobalSignal']))

        if self.data_controls['aCompCor']['enabled'].GetValue():
            selector['aCompCor'] = {}
            selector['aCompCor'].update(self.compile_selector_summary(self.data_controls['aCompCor']))
            selector['aCompCor']['summary']['method'] = 'DetrendPC'
            selector['aCompCor'].update(self.compile_selector_tissue_parameters(self.data_controls['aCompCor']))
            selector['aCompCor'].update(self.compile_selector_derivatives(self.data_controls['aCompCor']))

            selector['aCompCor']['tissues'] = []
            for t in range(self.data_controls['aCompCor']['tissues'].GetItemCount()):
                item = self.data_controls['aCompCor']['tissues'].GetItem(t).GetText()
                selector['aCompCor']['tissues'].append(selector_renaming_inverse[item])
            
        if self.data_controls['tCompCor']['enabled'].GetValue():
            selector['tCompCor'] = {}
            selector['tCompCor'].update(self.compile_selector_summary(self.data_controls['tCompCor']))
            selector['tCompCor'].update(self.compile_selector_derivatives(self.data_controls['tCompCor']))
            selector['tCompCor']['by_slice'] = self.data_controls['tCompCor']['by_slice'].GetValue()
            selector['tCompCor']['threshold'] = self.data_controls['tCompCor']['threshold'].GetValue()
            
        if self.data_controls['PolyOrt']['enabled'].GetValue():
            selector['PolyOrt'] = {
                'degree': int(self.data_controls['PolyOrt']['degree'].GetValue()),
            }
        
        if self.data_controls['Bandpass']['enabled'].GetValue():
            selector['Bandpass'] = {}

            bottom_frequency = self.data_controls['Bandpass']['bottom_frequency'].GetValue()
            if bottom_frequency:
                selector['Bandpass']['bottom_frequency'] = float(self.data_controls['Bandpass']['bottom_frequency'].GetValue())

            top_frequency = self.data_controls['Bandpass']['top_frequency'].GetValue()
            if top_frequency:
                selector['Bandpass']['top_frequency'] = float(self.data_controls['Bandpass']['top_frequency'].GetValue())

        if self.data_controls['Censor']['enabled'].GetValue():
            selector['Censor'] = {
                'method': str(self.data_controls['Censor']['method'].GetValue()),
                'number_of_previous_trs_to_censor': int(self.data_controls['Censor']['number_of_previous_trs_to_censor'].GetValue()) or 0,
                'number_of_subsequent_trs_to_censor': int(self.data_controls['Censor']['number_of_subsequent_trs_to_censor'].GetValue()) or 0,
                'thresholds': [],
            }

            if self.data_controls['Censor']['threshold']['fdj']['enabled'].GetValue():
                selector['Censor']['thresholds'].append({
                    'type': 'FD_J',
                    'value': self.parse_threshold(self.data_controls['Censor']['threshold']['fdj']['value']),
                })
            if self.data_controls['Censor']['threshold']['fdp']['enabled'].GetValue():
                selector['Censor']['thresholds'].append({
                    'type': 'FD_P',
                    'value': self.parse_threshold(self.data_controls['Censor']['threshold']['fdp']['value']),
                })
            if self.data_controls['Censor']['threshold']['dvars']['enabled'].GetValue():
                selector['Censor']['thresholds'].append({
                    'type': 'DVARS',
                    'value': self.parse_threshold(self.data_controls['Censor']['threshold']['dvars']['value']),
                })

        return selector

    def ok_click(self, event):
        evt = EditorOkEvent()
        wx.PostEvent(self, evt)
        self.Close()

    def cancel_click(self, event):
        evt = EditorCancelEvent()
        wx.PostEvent(self, evt)
        self.Close()


class NuisanceRegressionRegressorsGrid(wx.Panel):

    def __init__(self, parent, id=wx.ID_ANY, value=None, size=wx.DefaultSize):

        wx.Panel.__init__(self, parent, id=id, size=size)

        self.regressor_selectors = []

        sizer = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(sizer)

        self.scroll = scrolled.ScrolledPanel(self, size=(500, 250), style=wx.VSCROLL | wx.BORDER_SUNKEN)
        self.scroll.SetBackgroundColour(wx.WHITE)
        self.scroll.SetupScrolling(scroll_x=False, scroll_y=True)
        self.scroll_sizer = wx.GridBagSizer(vgap=5, hgap=5)
        self.scroll.SetSizer(self.scroll_sizer)
        sizer.Add(self.scroll, 1, wx.EXPAND | wx.ALL)

        bmp = wx.ArtProvider.GetBitmap(wx.ART_PLUS, wx.ART_OTHER, (16, 16))
        button = wx.BitmapButton(self, wx.ID_ANY, bmp)
        button.Bind(wx.EVT_BUTTON, self.add_regressor)
        sizer.Add(button, 0, wx.ALIGN_RIGHT | wx.ALL, border=5)

        self.render()

        self.Bind(wx.EVT_SHOW, self.render)
        self.Bind(wx.EVT_SIZE, self.render)

    def update(self, regressor_selectors):
        self.regressor_selectors = regressor_selectors
        self.render()

    def get_value(self):
        return self.regressor_selectors

    def add_regressor(self, event):
        editor = NuisanceRegressionRegressorEditor(self)

        def save_regressor(event):
            sel = editor.compile_selector()
            if sel:
                self.regressor_selectors.append(sel)
                self.render()

        editor.Bind(EVT_EDITOR_OK, save_regressor)
        editor.Show()

    def duplicate_regressor(self, event, regressor_i):

        cp = copy.deepcopy(self.regressor_selectors[regressor_i])
        regressor_i += 1 

        self.regressor_selectors = \
            self.regressor_selectors[0:regressor_i] + [cp] + self.regressor_selectors[regressor_i:]

        self.render()

        editor = NuisanceRegressionRegressorEditor(self, regressor_i)

        def save_regressor(event):
            self.regressor_selectors[regressor_i] = editor.compile_selector()
            self.render()

        editor.Bind(EVT_EDITOR_OK, save_regressor)
        editor.Show()

    def edit_regressor(self, event, regressor_i):
        editor = NuisanceRegressionRegressorEditor(self, regressor_i)

        def save_regressor(event):
            sel = editor.compile_selector()
            if sel:
                self.regressor_selectors[regressor_i] = sel
                self.render()

        editor.Bind(EVT_EDITOR_OK, save_regressor)
        editor.Show()

    def remove_regressor(self, event, regressor_i):
        del self.regressor_selectors[regressor_i]
        self.render()

    def render(self, event=None):

        for p in self.scroll.GetChildren():
            p.Destroy()
        self.scroll_sizer.SetRows(0)

        font = wx.Font(12, wx.DEFAULT, wx.NORMAL, wx.NORMAL)

        width, _ = self.scroll.GetSize()
        if event:
            width, _ = event.GetSize()

        for i, selectors in enumerate(self.regressor_selectors):

            row = self.scroll_sizer.Rows

            buttons_panel = wx.Panel(self.scroll)
            buttons_panel.SetBackgroundColour(wx.WHITE)
            buttons_sizer = wx.BoxSizer(wx.VERTICAL)
            buttons_panel.SetSizer(buttons_sizer)

            bmp = wx.ArtProvider.GetBitmap(wx.ART_EXECUTABLE_FILE, wx.ART_OTHER, (16, 16))
            button = wx.BitmapButton(buttons_panel, wx.ID_ANY, bmp)
            button.SetToolTip(wx.ToolTip("Edit regressors"))
            button.Bind(wx.EVT_BUTTON, (lambda i: lambda event: self.edit_regressor(event, i))(i))
            buttons_sizer.Add(button, 1, wx.CENTER | wx.ALL, border=3)

            bmp = wx.ArtProvider.GetBitmap(wx.ART_COPY, wx.ART_OTHER, (16, 16))
            button = wx.BitmapButton(buttons_panel, wx.ID_ANY, bmp)
            button.SetToolTip(wx.ToolTip("Duplicate regressors"))
            button.Bind(wx.EVT_BUTTON, (lambda i: lambda event: self.duplicate_regressor(event, i))(i))
            buttons_sizer.Add(button, 1, wx.CENTER | wx.ALL, border=3)
            
            bmp = wx.ArtProvider.GetBitmap(wx.ART_DELETE, wx.ART_OTHER, (16, 16))
            button = wx.BitmapButton(buttons_panel, wx.ID_ANY, bmp)
            button.SetToolTip(wx.ToolTip("Remove regressors"))
            button.Bind(wx.EVT_BUTTON, (lambda i: lambda event: self.remove_regressor(event, i))(i))
            buttons_sizer.Add(button, 1, wx.CENTER | wx.ALL, border=3)

            self.scroll_sizer.Add(buttons_panel, pos=(row, 0), flag=wx.ALIGN_CENTER_VERTICAL | wx.ALL, border=5)

            selectors_description = selectors_repr(selectors)
            description = wx.StaticText(self.scroll, wx.ID_ANY, selectors_description)
            description.Wrap(width - 200)
            description.SetFont(font)

            self.scroll_sizer.Add(description, pos=(row, 1), flag=wx.ALIGN_CENTER_VERTICAL | wx.ALL, border=5)

        self.scroll.Layout()
        self.Layout()


class NuisanceRegressionRegressors(Control):

    def __init__(self, parent):

        self.name = 'Regressors'
        self.default_values = []
        self.datatype = dtype.LDICT
        self.help = "Select which nuisance signal corrections to apply"

        self.type = self.__class__.__name__

        self.ctrl = NuisanceRegressionRegressorsGrid(
            parent,
            value={}
        )

    def set_value(self, value):
        self.ctrl.update(value)

    def get_validation(self):
        return False

    def get_selection(self):
        return self.ctrl.get_value()


class NuisanceRegression(wx.ScrolledWindow):

    def __init__(self, parent, counter = 0):
        wx.ScrolledWindow.__init__(self, parent)

        self.counter = counter
        self.page = GenericClass(self, "Nuisance Signal Regression Options")

        self.page.add(label="Run Nuisance Regression ",
                 control=control.CHOICE_BOX,
                 name='runNuisance',
                 type=dtype.LSTR,
                 comment="Run Nuisance Signal Regression",
                 values=["Off", "On", "On/Off"],
                 wkf_switch = True)

        self.page.add(label="Lateral Ventricles Mask\n(Standard Space) ",
                     control=control.COMBO_BOX,
                     name='lateral_ventricles_mask',
                     type=dtype.STR,
                     values="$FSLDIR/data/atlases/HarvardOxford/HarvardOxford-lateral-ventricles-thr25-2mm.nii.gz",
                     comment="Standard Lateral Ventricles Binary Mask")

        self.page.add(label="Select Regressors\nand Censoring:",
                      control=NuisanceRegressionRegressors(self),
                      name="Regressors",
                      comment="List of regressors, censoring and detrending"
                              "that you want to run",
                      type=dtype.LDICT)

        self.page.set_sizer()
        self.page.flexSizer.AddGrowableRow(2)
        self.page.flexSizer.Layout()
        self.GetSizer().Layout()

        if hasattr(parent, 'get_page_list'):
            parent.get_page_list().append(self)

    def get_counter(self):
            return self.counter


class MedianAngleCorrection(wx.ScrolledWindow):

    def __init__(self, parent, counter = 0):
        wx.ScrolledWindow.__init__(self, parent)

        self.counter = counter

        self.page = GenericClass(self, "Median Angle Correction Options")

        self.page.add(label="Run Median Angle Correction ",
                 control=control.CHOICE_BOX,
                 name='runMedianAngleCorrection',
                 type=dtype.LSTR,
                 comment="Correct for the global signal using Median Angle Correction.",
                 values=["Off","On","On/Off"],
                 wkf_switch = True)

        self.page.add(label= "Target Angle (degrees) ",
                      control = control.TEXT_BOX,
                      name = "targetAngleDeg",
                      type = dtype.LNUM,
                      values = "90",
                      validator = CharValidator("no-alpha"),
                      comment = "Target angle used during Median Angle Correction.")

        self.page.set_sizer()
        parent.get_page_list().append(self)

    def get_counter(self):
            return self.counter
