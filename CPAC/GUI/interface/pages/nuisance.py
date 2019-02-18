# -*- coding: utf-8 -*-

import copy
import wx
import wx.html
import wx.lib.newevent
from wx.lib.mixins.listctrl import CheckListCtrlMixin, ListCtrlAutoWidthMixin
from ..utils.generic_class import GenericClass, Control
from ..utils.constants import control, dtype
from ..utils.validator import CharValidator
import pkg_resources as p
from collections import defaultdict


regressor_selectors = [
    {
        'CerebrospinalFluid': {'erode_mask': True,
                               'extraction_resolution': 2,
                               'summary': {'components': 5, 'method': 'PC'}},
        'GlobalSignal': {'include_delayed': True,
                         'include_delayed_squared': True,
                         'include_squared': True,
                         'summary': 'Mean'},
        'Motion': {'include_delayed': True,
                   'include_delayed_squared': True,
                   'include_squared': True},
        'WhiteMatter': {'extraction_resolution': 2,
                        'summary': {'components': 5, 'method': 'PC'}},
        'aCompCor': {'extraction_resolution': 2,
                     'summary': {'components': 5, 'method': 'PC'},
                     'tissues': ['WhiteMatter', 'CerebrospinalFluid']},
        'tCompCor': {'by_slice': True,
                     'summary': {'components': 5, 'method': 'PC'},
                     'threshold': '1.5SD'},
        'PolyOrt': {'degree': 2},
        'Bandpass': {
            'bottom_frequency': 0.01,
            'top_frequency': 0.1,
        },
        'Censor': {'method': 'Interpolate',
                   'thresholds': [
                       # {'type': 'FD', 'value': 0.5},
                       {'type': 'DVARS', 'value': 17}
                   ]},
    },
    {
        'Motion': {'include_delayed': True,
                   'include_delayed_squared': True,
                   'include_squared': True},
        'Censor': {'method': 'Interpolate',
                   'thresholds': [{'type': 'FD', 'value': 0.5},
                                  {'type': 'DVARS', 'value': 17}]},
    },
]

def find(lst, lmbd):
    return next((t for t in lst if lmbd(t)), None)

def findi(lst, lmbd):
    return next((i for i, t in enumerate(lst) if lmbd(t)), None)

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
        regressor = selectors[reg]

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
                if regressor['summary']['method'] == 'PC':
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

        representation += "{} ({})\n".format(method, ", ".join(threshs))

    bandpass = selectors.get('Bandpass')
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

    polynomial = selectors.get('PolyOrt')
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
    'PolyOrt': 'Poly Regression',
    'Bandpass': 'Bandpass',
    'Censor': 'Censor',
}

selector_renaming_inverse = \
    dict(zip(*zip(*selector_renaming.items())[::-1]))

selector_summary_method_renaming = {
    'PC': 'Principal Components',
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

    def __init__(self, parent, selectors_index):
        wx.Frame.__init__(self, parent,
                          title="Edit Selectors",
                          size=(650, 400))

        self.selectors_index = selectors_index
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

        self.editor = wx.ScrolledWindow(editor, style=wx.VSCROLL | wx.BORDER_SUNKEN)
        self.editor.SetScrollRate(0, 5)
        self.editor.EnableScrolling(False, True)
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
        button_ok.Bind(wx.EVT_BUTTON, self.onButtonClick)
        button_sizer.Add(button_ok, 0, wx.ALIGN_CENTER)

        button_cancel = wx.Button(button_panel, -1, 'Cancel', size=(90, 30))
        button_cancel.Bind(wx.EVT_BUTTON, self.onButtonClick)
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

        control = wx.ComboBox(summary, choices=selector_summary_method_renaming.values())
        control.SetValue(selector_summary_method_renaming[selector_summary_method])
        summary_sizer.Add(control, flag=wx.EXPAND | wx.ALL)

        self.add_to_new_row(label, summary)
        controls['summary']['method'] = control

        def control_update(event=None):
            if selector_summary_method_renaming_inverse.get(control.GetValue()) == 'PC':
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
        self.data_controls['aCompCor']['components'] = acompcor_components_control
        
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
        self.data_controls['tCompCor']['components'] = tcompcor_components_control

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


        self.titles += [self.render_selectors_title('Detrending', main=True)]

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
        censoring_method_control = wx.ComboBox(self.editor, choices=selector_censor_method_renaming.values())
        if censoring_method:
            censoring_method_control.SetValue(selector_censor_method_renaming[censoring_method])
        self.add_to_new_row(censoring_method_label, censoring_method_control)
        self.data_controls['Censor']['method'] = censoring_method_control

        # Framewise displacement control
        censoring_thresholds_fd = find(censoring_thresholds, lambda t:  t.get('type') == 'FD') or {}
        has_censoring_thresholds_fd = bool(censoring_thresholds_fd)

        censoring_thresholds_fd_label = wx.StaticText(self.editor, label='Framewise Displacement')
        censoring_thresholds_fd_control = wx.CheckBox(self.editor)
        censoring_thresholds_fd_control.SetValue(has_censoring_thresholds_fd)
        self.add_to_new_row(censoring_thresholds_fd_label, censoring_thresholds_fd_control)
        self.data_controls['Censor']['threshold']['fd']['enabled'] = censoring_thresholds_fd_control

        censoring_thresholds_fd_value_label = wx.StaticText(self.editor, label='Threshold')
        censoring_thresholds_fd_value_control = wx.TextCtrl(self.editor)
        censoring_thresholds_fd_value_control.SetValue(str(censoring_thresholds_fd.get('value', 0.0)))
        self.add_to_new_row(censoring_thresholds_fd_value_label, censoring_thresholds_fd_value_control)
        self.data_controls['Censor']['threshold']['fd']['value'] = censoring_thresholds_fd_value_control

        def censoring_thresholds_fd_update(event=None):
            if censoring_thresholds_fd_control.GetValue():
                ctrl_enable(censoring_thresholds_fd_value_control)
            else:
                ctrl_disable(censoring_thresholds_fd_value_control)

        censoring_thresholds_fd_control.Bind(wx.EVT_CHECKBOX, censoring_thresholds_fd_update)
        censoring_thresholds_fd_control.Bind(EVT_ENABLE, censoring_thresholds_fd_update)
        censoring_thresholds_fd_update()

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

        def censoring_thresholds_dvars_active(event=None):
             print("Active!!!")

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
        ]

        for k in valid_regressors:
            self.tree.AppendItem(root_regressors, selector_renaming[k])

        valid_detrending = [
            'PolyOrt',
            'Bandpass',
        ]

        for k in valid_detrending:
            if k in self.selectors:
                self.tree.AppendItem(root_detrending, selector_renaming[k])

        self.tree.ExpandAll()

    def onButtonClick(self,event):
        self.Close()


class NuisanceRegressionRegressorsGrid(wx.Panel):

    def __init__(self, parent, id=wx.ID_ANY, value=None, size=wx.DefaultSize):

        wx.Panel.__init__(self, parent, id=id)

        self.regressor_selectors = regressor_selectors

        self.scroll = wx.ScrolledWindow(self, size=(500, 250), style=wx.VSCROLL | wx.BORDER_SUNKEN)
        self.scroll.SetScrollRate(0, 10)
        self.scroll.EnableScrolling(False, True)
        self.scroll.SetBackgroundColour(wx.WHITE)

        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.scroll_sizer = wx.BoxSizer(wx.VERTICAL)

        self.render()

        self.scroll.SetSizer(self.scroll_sizer)
        self.sizer.Add(self.scroll, 0, wx.EXPAND | wx.ALL)
        self.SetSizer(self.sizer)

        self.Bind(wx.EVT_SHOW, self.render)
        self.Bind(wx.EVT_SIZE, self.render)

    def edit_regressor(self, event, regressor_i):
        editor = NuisanceRegressionRegressorEditor(self, regressor_i)
        editor.Show()

    def remove_regressor(self, event, regressor_i):
        del self.regressor_selectors[regressor_i]
        self.render()

    def render(self, event=None):

        for p in self.scroll.GetChildren():
            p.Destroy()

        font = wx.Font(12, wx.DEFAULT, wx.NORMAL, wx.NORMAL)

        width, _ = self.scroll.GetSize()
        if event:
            width, _ = event.GetSize()

        for i, selectors in enumerate(self.regressor_selectors):

            item_panel = wx.Panel(self.scroll, style=wx.BORDER_RAISED)
            item_panel.SetBackgroundColour(wx.WHITE)
            item_sizer = wx.BoxSizer(wx.HORIZONTAL)
            item_panel.SetSizer(item_sizer)

            buttons_panel = wx.Panel(item_panel)
            buttons_panel.SetBackgroundColour(wx.WHITE)
            buttons_sizer = wx.BoxSizer(wx.VERTICAL)
            buttons_panel.SetSizer(buttons_sizer)

            bmp = wx.ArtProvider.GetBitmap(wx.ART_COPY, wx.ART_OTHER, (16, 16))
            button = wx.BitmapButton(buttons_panel, wx.ID_ANY, bmp)
            button.Bind(wx.EVT_BUTTON, (lambda i: lambda event: self.edit_regressor(event, i))(i))
            buttons_sizer.Add(button, 1, wx.CENTER | wx.ALL, border=5)

            bmp = wx.ArtProvider.GetBitmap(wx.ART_DELETE, wx.ART_OTHER, (16, 16))
            button = wx.BitmapButton(buttons_panel, wx.ID_ANY, bmp)
            button.Bind(wx.EVT_BUTTON, (lambda i: lambda event: self.remove_regressor(event, i))(i))
            buttons_sizer.Add(button, 1, wx.CENTER | wx.ALL, border=5)

            item_sizer.Add(buttons_panel, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, border=5)

            selectors_description = selectors_repr(selectors)
            description = wx.StaticText(item_panel, wx.ID_ANY, selectors_description)
            description.Wrap(width - 200)
            description.SetFont(font)
    
            item_sizer.Add(description, 1, wx.ALIGN_CENTER_VERTICAL | wx.ALL | wx.ST_ELLIPSIZE_END, border=5)

            self.scroll_sizer.Add(item_panel, 0, wx.EXPAND | wx.ALL, 2)

        self.scroll.Layout()
        self.Layout()


class NuisanceRegressionRegressors(Control):

    def __init__(self, parent):

        self.name = 'Regressors'
        self.default_values = []

        self.ctrl = NuisanceRegressionRegressorsGrid(
            parent,
            value=regressor_selectors
        )


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


def test_nuisance():

    app = wx.App(False)
    frame = wx.Frame(None, wx.ID_ANY, "Nuisance Regression")
    
    # editor = NuisanceRegression(frame)
    # frame.Show(True)
    # frame.Bind(wx.EVT_CLOSE, lambda event: app.Destroy())

    frame.regressor_selectors = regressor_selectors
    editor = NuisanceRegressionRegressorEditor(frame, 0)
    editor.Show(True)
    editor.Bind(wx.EVT_CLOSE, lambda event: app.Destroy())

    app.MainLoop()

