# -*- coding: utf-8 -*-

import copy
import wx
import wx.html
from ..utils.generic_class import GenericClass, Control
from ..utils.constants import control, dtype
from ..utils.validator import CharValidator
import pkg_resources as p


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
                        {'type': 'FD', 'value': 0.5},
                        {'type': 'DVARS', 'value': 17}
                    ]},
    },
    {
        'CerebrospinalFluid': {'erode_mask': True,
                               'extraction_resolution': 2,
                               'summary': {'components': 5, 'method': 'PC'}},
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
                   'thresholds': [{'type': 'FD', 'value': 0.5},
                                  {'type': 'DVARS', 'value': 17}]},
    },
]

def find(lst, lmbd):
    return next((t for t in lst if lmbd(t)), None)

def findi(lst, lmbd):
    return next((i for i, t in enumerate(lst) if lmbd(t)), None)

def ctrl_enable(ctrl):
    ctrl.SetEditable(True)
    ctrl.SetBackgroundColour(wx.WHITE)

def ctrl_disable(ctrl):
    color = wx.SystemSettings.GetColour(wx.SYS_COLOUR_BACKGROUND)
    ctrl.SetEditable(False)
    ctrl.SetBackgroundColour(color)

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
        if regressor.get('include_delayed'):
            terms += ["{}{}".format(renamed, 'ₜ₋₁')]
        if regressor.get('include_squared'):
            terms += ["{}{}".format(renamed, '²')]
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


class NuisanceRegressionRegressorEditor(wx.Frame):

    def __init__(self, parent, selectors_index):
        wx.Frame.__init__(self, parent,
                          title="Edit regressors",
                          size=(520,90))

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
                                style=wx.TR_DEFAULT_STYLE |
                                      wx.TR_FULL_ROW_HIGHLIGHT |
                                      wx.TR_LINES_AT_ROOT |
                                      wx.TR_HIDE_ROOT)
        self.Bind(wx.EVT_TREE_SEL_CHANGED, self.scroll_selector, self.tree)

        # Render Editor
        editor = wx.Panel(panel)
        editor_sizer = wx.BoxSizer(wx.VERTICAL)
        editor.SetSizer(editor_sizer)

        self.editor = wx.ScrolledWindow(editor, style=wx.VSCROLL | wx.BORDER_SUNKEN)
        self.editor.SetScrollRate(0, 5)
        self.editor.EnableScrolling(False, True)
        self.editor.SetBackgroundColour(wx.WHITE)
        self.editor.SetSizer(wx.BoxSizer(wx.VERTICAL))

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

    def render_selectors_title(self, title):
        sizer = self.editor.GetSizer()
        font = wx.Font(12, wx.DEFAULT, wx.NORMAL, wx.FONTWEIGHT_BOLD)
        description = wx.StaticText(self.editor, wx.ID_ANY, title)
        description.SetFont(font)
        sizer.Add(description, flag=wx.EXPAND | wx.ALL, border=5)
        return description

    def render_selectors(self):
        
        parts = [
            'Motion',
            'CerebrospinalFluid',
            'WhiteMatter',
            'GreyMatter',
            'GlobalSignal',
            'aCompCor',
            'tCompCor',
            'PolyOrt',
            'Bandpass',
            'Censor',
        ]

        sizer = self.editor.GetSizer()

        self.titles = []

        self.selectors = {
            'PolyOrt': {
                'degree': 2
            },
            'Bandpass': {
                'bottom_frequency': 0.01,
                'top_frequency': 0.1,
            },
            'Censor': {'method': 'Interpolate',
                       'thresholds': [{'type': 'FD', 'value': 0.5},
                                      {'type': 'DVARS', 'value': 17}]},
        }

        self.titles += [self.render_selectors_title('Detrending')]


        self.titles += [self.render_selectors_title('Poly Regression')]
        
        has_poly = 'PolyOrt' in self.selectors

        poly_panel = wx.Panel(self.editor)
        poly_panel.SetBackgroundColour(wx.WHITE)
        poly_panel_sizer = wx.FlexGridSizer(cols=2)
        poly_panel_sizer.AddGrowableCol(1)
        poly_panel.SetSizer(poly_panel_sizer)
        sizer.Add(poly_panel, flag=wx.EXPAND | wx.ALL, border=5)

        poly_enabled_label = wx.StaticText(poly_panel, wx.ID_ANY, label='Enabled', style=wx.ALIGN_RIGHT)
        poly_enabled_control = wx.CheckBox(poly_panel, wx.ID_ANY)
        poly_enabled_control.SetValue(has_poly)
        poly_panel_sizer.Add(poly_enabled_label, 1, wx.EXPAND | wx.ALL, border=5)
        poly_panel_sizer.Add(poly_enabled_control, 2, wx.EXPAND | wx.ALL, border=5)
        

        poly_degree_label = wx.StaticText(poly_panel, wx.ID_ANY, label='Degree', style=wx.ALIGN_RIGHT)
        poly_degree_control = wx.TextCtrl(poly_panel, id=wx.ID_ANY, value='0')
        poly_panel_sizer.Add(poly_degree_label, 1, wx.EXPAND | wx.ALL, border=5)
        poly_panel_sizer.Add(poly_degree_control, 2, wx.EXPAND | wx.ALL, border=5)
        if has_poly:
            poly_degree_control.SetValue(str(self.selectors['PolyOrt']['degree']))
        else:
            ctrl_disable(poly_degree_control)

        def poly_enabled_update(event):
            if poly_enabled_control.GetValue():
                ctrl_enable(poly_degree_control)
                self.selectors['PolyOrt'] = {
                    'degree': 2,
                }
                poly_degree_control.SetValue(str(self.selectors['PolyOrt']['degree']))
            else:
                ctrl_disable(poly_degree_control)
                del self.selectors['PolyOrt']

        poly_enabled_control.Bind(
            wx.EVT_CHECKBOX, poly_enabled_update
        )


        self.titles += [self.render_selectors_title('Bandpass')]

        has_bandpass = 'Bandpass' in self.selectors

        bandpass_freq_panel = wx.Panel(self.editor)
        bandpass_freq_panel.SetBackgroundColour(wx.WHITE)
        bandpass_freq_panel_sizer = wx.FlexGridSizer(cols=2)
        bandpass_freq_panel_sizer.AddGrowableCol(1)
        bandpass_freq_panel.SetSizer(bandpass_freq_panel_sizer)
        sizer.Add(bandpass_freq_panel, flag=wx.EXPAND | wx.ALL, border=5)

        bandpass_freq_enabled_label = wx.StaticText(bandpass_freq_panel, wx.ID_ANY, label='Enabled', style=wx.ALIGN_RIGHT)
        bandpass_freq_enabled_control = wx.CheckBox(bandpass_freq_panel, wx.ID_ANY)
        bandpass_freq_enabled_control.SetValue(has_bandpass)
        bandpass_freq_panel_sizer.Add(bandpass_freq_enabled_label, 1, wx.EXPAND | wx.ALL, border=5)
        bandpass_freq_panel_sizer.Add(bandpass_freq_enabled_control, 2, wx.EXPAND | wx.ALL, border=5)

        bandpass_freq_bottom_label = wx.StaticText(bandpass_freq_panel, wx.ID_ANY, label='Bottom Frequency', style=wx.ALIGN_RIGHT)
        bandpass_freq_bottom_control = wx.TextCtrl(bandpass_freq_panel, id=wx.ID_ANY, value='0.0')
        bandpass_freq_panel_sizer.Add(bandpass_freq_bottom_label, 1, wx.EXPAND | wx.ALL, border=5)
        bandpass_freq_panel_sizer.Add(bandpass_freq_bottom_control, 2, wx.EXPAND | wx.ALL, border=5)
        if has_bandpass:
            bandpass_freq_bottom_control.SetValue(str(self.selectors['Bandpass']['bottom_frequency']))
        else:
            ctrl_disable(bandpass_freq_bottom_control)
        

        bandpass_freq_top_label = wx.StaticText(bandpass_freq_panel, wx.ID_ANY, label='Top Frequency', style=wx.ALIGN_RIGHT)
        bandpass_freq_top_control = wx.TextCtrl(bandpass_freq_panel, id=wx.ID_ANY, value='0.0')
        bandpass_freq_panel_sizer.Add(bandpass_freq_top_label, 1, wx.EXPAND | wx.ALL, border=5)
        bandpass_freq_panel_sizer.Add(bandpass_freq_top_control, 2, wx.EXPAND | wx.ALL, border=5)
        if has_bandpass:
            bandpass_freq_top_control.SetValue(str(self.selectors['Bandpass']['top_frequency']))
        else:
            ctrl_disable(bandpass_freq_top_control)
        

        def bandpass_freq_enabled_update(event):
            if bandpass_freq_enabled_control.GetValue():
                ctrl_enable(bandpass_freq_bottom_control)
                ctrl_enable(bandpass_freq_top_control)
                self.selectors['Bandpass'] = {
                    'bottom_frequency': 0.01,
                    'top_frequency': 0.1,
                }
                bandpass_freq_bottom_control.SetValue(str(self.selectors['Bandpass']['bottom_frequency']))
                bandpass_freq_top_control.SetValue(str(self.selectors['Bandpass']['top_frequency']))
            else:
                ctrl_disable(bandpass_freq_bottom_control)
                ctrl_disable(bandpass_freq_top_control)
                del self.selectors['Bandpass']

        bandpass_freq_enabled_control.Bind(
            wx.EVT_CHECKBOX, bandpass_freq_enabled_update
        )


        self.titles += [self.render_selectors_title('Censoring')]

        censoring = self.selectors.get('Censor', {})
        censoring_thresholds = censoring.get('thresholds', [])

        # Framewise displacement control
        censoring_thresholds_fd = find(censoring_thresholds, lambda t:  t.get('type') == 'FD')
        has_censoring_thresholds_fd = censoring_thresholds_fd is not None
        censoring_thresholds_fd_checkbox = wx.CheckBox(self.editor, wx.ID_ANY, 'Framewise Displacement')
        censoring_thresholds_fd_checkbox.SetValue(has_censoring_thresholds_fd)
        sizer.Add(censoring_thresholds_fd_checkbox, flag=wx.EXPAND | wx.ALL, border=5)

        censoring_thresholds_fd_thresh_panel = wx.Panel(self.editor)
        censoring_thresholds_fd_thresh_panel.SetBackgroundColour(wx.WHITE)
        censoring_thresholds_fd_thresh_panel_sizer = wx.BoxSizer(wx.HORIZONTAL)
        censoring_thresholds_fd_thresh_panel.SetSizer(censoring_thresholds_fd_thresh_panel_sizer)
        sizer.Add(censoring_thresholds_fd_thresh_panel, flag=wx.EXPAND | wx.ALL, border=5)

        censoring_thresholds_fd_thresh_label = wx.StaticText(censoring_thresholds_fd_thresh_panel, wx.ID_ANY, 'Threshold')
        censoring_thresholds_fd_thresh_panel_sizer.Add(censoring_thresholds_fd_thresh_label, 1, wx.CENTER | wx.ALL, border=5)

        censoring_thresholds_fd_thresh_control = wx.TextCtrl(censoring_thresholds_fd_thresh_panel, id=wx.ID_ANY, value=str(censoring_thresholds_fd['value']))
        censoring_thresholds_fd_thresh_panel_sizer.Add(censoring_thresholds_fd_thresh_control, 2, wx.CENTER | wx.ALL, border=5)
        if not has_censoring_thresholds_fd:
            ctrl_disable(censoring_thresholds_fd_thresh_control)

        def censoring_thresholds_fd_update(event):
            if 'Censor' not in self.selectors:
                self.selectors['Censor'] = {}
            if 'thresholds' not in self.selectors['Censor']:
                self.selectors['Censor']['thresholds'] = []

            fd = findi(self.selectors['Censor']['thresholds'], lambda t:  t.get('type') == 'FD')
            if censoring_thresholds_fd_checkbox.GetValue():
                if fd is None:
                    self.selectors['Censor']['thresholds'].append({"type": "FD", "value": 0.0})
                    ctrl_enable(censoring_thresholds_fd_thresh_control)
            else:
                if fd is not None:
                    del self.selectors['Censor']['thresholds'][fd]
                ctrl_disable(censoring_thresholds_fd_thresh_control)

        censoring_thresholds_fd_checkbox.Bind(
            wx.EVT_CHECKBOX, censoring_thresholds_fd_update
        )

        # DVARS control
        censoring_thresholds_dvars = find(censoring_thresholds, lambda t:  t.get('type') == 'DVARS')
        has_censoring_thresholds_dvars = censoring_thresholds_dvars is not None
        censoring_thresholds_dvars_checkbox = wx.CheckBox(self.editor, wx.ID_ANY, 'DVARS')
        censoring_thresholds_dvars_checkbox.SetValue(has_censoring_thresholds_dvars)
        sizer.Add(censoring_thresholds_dvars_checkbox, flag=wx.EXPAND | wx.ALL, border=5)

        censoring_thresholds_dvars_thresh_panel = wx.Panel(self.editor)
        censoring_thresholds_dvars_thresh_panel.SetBackgroundColour(wx.WHITE)
        censoring_thresholds_dvars_thresh_panel_sizer = wx.BoxSizer(wx.HORIZONTAL)
        censoring_thresholds_dvars_thresh_panel.SetSizer(censoring_thresholds_dvars_thresh_panel_sizer)
        sizer.Add(censoring_thresholds_dvars_thresh_panel, flag=wx.EXPAND | wx.ALL, border=5)

        censoring_thresholds_dvars_thresh_label = wx.StaticText(censoring_thresholds_dvars_thresh_panel, wx.ID_ANY, 'Threshold')
        censoring_thresholds_dvars_thresh_panel_sizer.Add(censoring_thresholds_dvars_thresh_label, 1, wx.CENTER | wx.ALL, border=5)

        censoring_thresholds_dvars_thresh_control = wx.TextCtrl(censoring_thresholds_dvars_thresh_panel, id=wx.ID_ANY, value=str(censoring_thresholds_dvars['value']))
        censoring_thresholds_dvars_thresh_panel_sizer.Add(censoring_thresholds_dvars_thresh_control, 2, wx.CENTER | wx.ALL, border=5)
        if not has_censoring_thresholds_dvars:
            ctrl_disable(censoring_thresholds_dvars_thresh_control)

        def censoring_thresholds_dvars_update(event):
            if 'Censor' not in self.selectors:
                self.selectors['Censor'] = {}
            if 'thresholds' not in self.selectors['Censor']:
                self.selectors['Censor']['thresholds'] = []

            dvars = findi(self.selectors['Censor']['thresholds'], lambda t:  t.get('type') == 'DVARS')
            if censoring_thresholds_dvars_checkbox.GetValue():
                if dvars is None:
                    self.selectors['Censor']['thresholds'].append({"type": "DVARS", "value": 0.0})
                    ctrl_enable(censoring_thresholds_dvars_thresh_control)
            else:
                if dvars is not None:
                    del self.selectors['Censor']['thresholds'][dvars]
                ctrl_disable(censoring_thresholds_dvars_thresh_control)

        censoring_thresholds_dvars_checkbox.Bind(
            wx.EVT_CHECKBOX, censoring_thresholds_dvars_update
        )



    def scroll_selector(self, event):
        item_name = str(self.tree.GetItemText(event.GetItem()))
        for title in self.titles:
            if item_name == title.GetLabelText():
                _, sppy = self.editor.GetScrollPixelsPerUnit()
                sy = self.editor.GetScrollPos(wx.VERTICAL) * sppy
                _, y = title.GetPosition().Get()
                self.editor.Scroll(0, int((sy + y) / sppy) - 5)

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
            if k in self.selectors:
                self.tree.AppendItem(root_regressors, selector_renaming[k])

        valid_detrending = [
            'PolyOrt',
            'Bandpass',
        ]

        for k in valid_detrending:
            if k in self.selectors:
                self.tree.AppendItem(root_detrending, selector_renaming[k])

    def onButtonClick(self,event):
        self.Close()


class NuisanceRegressionRegressorsGrid(wx.Panel):

    def __init__(self, parent, id, regressor_selectors, size):

        wx.Panel.__init__(self, parent, id=id)

        self.regressor_selectors = regressor_selectors

        self.scroll = wx.ScrolledWindow(self, size=(500, 250), style=wx.VSCROLL | wx.BORDER_SUNKEN)
        self.scroll.SetScrollRate(10,10)
        self.scroll.EnableScrolling(True,True)
        self.scroll.SetBackgroundColour(wx.WHITE)

        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.scroll_sizer = wx.BoxSizer(wx.VERTICAL)
        
        self.rows = {}
        self.render()

        self.scroll.SetSizer(self.scroll_sizer)
        self.sizer.Add(self.scroll, 0, wx.EXPAND | wx.ALL)
        self.SetSizer(self.sizer)

    def edit_regressor(self, event, regressor_i):
        editor = NuisanceRegressionRegressorEditor(self, regressor_i)
        editor.Show()

    def remove_regressor(self, event, regressor_i):
        del self.regressor_selectors[regressor_i]
        self.render()

    def render(self):

        # why not draw everything again
        for p in self.scroll.GetChildren():
            p.Destroy()
        # for i in self.rows:
        #     for el in self.rows[i]:
        #         el.Destroy()

        font = wx.Font(12, wx.DEFAULT, wx.NORMAL, wx.NORMAL)

        print(len(self.regressor_selectors))

        self.rows = {}

        for i, selectors in enumerate(self.regressor_selectors):

            self.rows[i] = []

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
            
            selectors_description = selectors_repr(selectors)
            description = wx.StaticText(item_panel, wx.ID_ANY, selectors_description)
            description.SetFont(font)

            item_sizer.Add(buttons_panel, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, border=5)
            item_sizer.Add(description, 1, wx.ALIGN_CENTER_VERTICAL | wx.ALL, border=5)

            self.scroll_sizer.Add(item_panel, 0, wx.EXPAND | wx.ALL, 2)

            # self.rows[i].append(item_panel)


        # self.scroll_sizer.Layout()
        # self.topSizer.Layout()

        self.scroll.Layout()
        self.Layout()


class NuisanceRegressionRegressors(Control):
    
    def __init__(self, parent):

        self.name = 'Regressors'
        self.default_values = []

        self.ctrl = NuisanceRegressionRegressorsGrid(
            parent, id=wx.ID_ANY,
            regressor_selectors=regressor_selectors,
            size=wx.DefaultSize
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
    # NuisanceRegression(frame)
    # frame.Show(True)

    frame.regressor_selectors = regressor_selectors
    editor = NuisanceRegressionRegressorEditor(frame, 0)
    editor.Show(True)
    editor.Bind(wx.EVT_CLOSE, lambda event: app.Destroy())
    
    app.MainLoop()

