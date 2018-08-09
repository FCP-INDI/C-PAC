import wx
import wx.html
from ..utils.generic_class import GenericClass
from ..utils.constants import control, dtype
from ..utils.custom_control import ListBoxComboEditor
import os
import pkg_resources as p


class TypeOption:

    def __init__(self, type, parameters):
        assert type in self.labels
        self.type = type
        self.parameters = parameters

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        params = ''
        if self.parameters:
            fparams = ['%s: %s' % (k, v) for k, v in self.parameters.items()]
            params = ' [%s]' %  ("; ".join(fparams))

        return self.labels[self.type] + params

    def __eq__(self, other): 
        return self.__dict__ == other.__dict__


class RoiOption(TypeOption):
    
    labels = {
        'dict_learning': 'Dictionary Learning',
        'k_means': 'k-Means',
        'ward': 'Ward',
        'group_ica': 'Group ICA',
    }

class ConnectivityOption(TypeOption):
    
    labels = {
        'correlation': 'Correlation',
        'partial_correlation': 'Partial Correlation',
        'tangent': 'Tangent',
    }

class ClassifierOption(TypeOption):
    
    labels = {
        'svm': 'Support Vector Machine',
        'ridge': 'Ridge',
        'logistic': 'Logistic Regresson',
        'random_forest': 'Random Forest',
        'knn': 'K Nearest Neighbors',
        'naive_bayes': 'Gaussian Naive Bayes',
    }


class RoiEditor(ListBoxComboEditor):

    def __init__(self, parent):
        super(RoiEditor, self).__init__(
            parent=parent, title="Region of Interest Selector",
            size=(300, 300)
        )
        
        self.parent = parent

        panel = wx.Panel(self)
        sizer = wx.BoxSizer(wx.VERTICAL)
        flexsizer = wx.FlexGridSizer(cols=2, hgap=10, vgap=15)

        labels = RoiOption.labels

        type_label = wx.StaticText(panel, -1, label='Type')
        self.type = wx.Choice(panel, wx.ID_ANY,
                              validator=wx.DefaultValidator,
                              choices=[labels[k] for k in sorted(labels.keys())])

        flexsizer.Add(type_label)
        flexsizer.Add(self.type)

        button = wx.Button(panel, -1, 'OK', size=(90, 30))
        button.Bind(wx.EVT_BUTTON, self.onFinish)
        sizer.Add(flexsizer, 1, wx.EXPAND | wx.ALL, 10)
        sizer.Add(button, 0, wx.ALIGN_CENTER)
        panel.SetSizer(sizer)

    def onFinish(self,event):
        parent = self.Parent
        if self.type.GetSelection() > -1:
            labels = sorted(RoiOption.labels.keys())
            val = labels[self.type.GetSelection()]

            option = RoiOption(val, {})
            
            parent.listbox.Append(str(option), option)
            # TODO check appended item
            # TODO check for duplicates

            self.Close()


class ConnectivityEditor(ListBoxComboEditor):

    def __init__(self, parent):
        super(ConnectivityEditor, self).__init__(
            parent=parent, title="Region of Interest Selector",
            size=(300, 300)
        )
        
        self.parent = parent

        panel = wx.Panel(self)
        sizer = wx.BoxSizer(wx.VERTICAL)
        flexsizer = wx.FlexGridSizer(cols=2, hgap=10, vgap=15)

        labels = ConnectivityOption.labels

        type_label = wx.StaticText(panel, -1, label='Type')
        self.type = wx.Choice(panel, wx.ID_ANY,
                              validator=wx.DefaultValidator,
                              choices=[labels[k] for k in sorted(labels.keys())])

        flexsizer.Add(type_label)
        flexsizer.Add(self.type)

        button = wx.Button(panel, -1, 'OK', size=(90, 30))
        button.Bind(wx.EVT_BUTTON, self.onFinish)
        sizer.Add(flexsizer, 1, wx.EXPAND | wx.ALL, 10)
        sizer.Add(button, 0, wx.ALIGN_CENTER)
        panel.SetSizer(sizer)

    def onFinish(self,event):
        parent = self.Parent
        if self.type.GetSelection() > -1:
            labels = sorted(ConnectivityOption.labels.keys())
            val = labels[self.type.GetSelection()]

            option = ConnectivityOption(val, {})
            
            parent.listbox.Append(str(option), option)
            # TODO check appended item
            # TODO check for duplicates

            self.Close()


class ClassifierEditor(ListBoxComboEditor):

    def __init__(self, parent):
        super(ClassifierEditor, self).__init__(
            parent=parent, title="Region of Interest Selector",
            size=(300, 300)
        )
        
        self.parent = parent

        panel = wx.Panel(self)
        sizer = wx.BoxSizer(wx.VERTICAL)
        flexsizer = wx.FlexGridSizer(cols=2, hgap=10, vgap=15)

        labels = ClassifierOption.labels

        type_label = wx.StaticText(panel, -1, label='Type')
        self.type = wx.Choice(panel, wx.ID_ANY,
                              validator=wx.DefaultValidator,
                              choices=[labels[k] for k in sorted(labels.keys())])

        flexsizer.Add(type_label)
        flexsizer.Add(self.type)

        button = wx.Button(panel, -1, 'OK', size=(90, 30))
        button.Bind(wx.EVT_BUTTON, self.onFinish)
        sizer.Add(flexsizer, 1, wx.EXPAND | wx.ALL, 10)
        sizer.Add(button, 0, wx.ALIGN_CENTER)
        panel.SetSizer(sizer)

    def onFinish(self,event):
        parent = self.Parent
        if self.type.GetSelection() > -1:
            labels = sorted(ClassifierOption.labels.keys())
            val = labels[self.type.GetSelection()]

            option = ClassifierOption(val, {})
            
            parent.listbox.Append(str(option), option)
            # TODO check appended item
            # TODO check for duplicates

            self.Close()

    
class ConnectomeSettings(wx.ScrolledWindow):
    def __init__(self, parent, counter = 0):
        wx.ScrolledWindow.__init__(self, parent)
        
        self.counter = counter
                
        self.page = GenericClass(self, "Connectome studies")

        self.page.add(label="Run Connectome?",
                      control=control.CHOICE_BOX,
                      name="connectome.run",
                      type=dtype.BOOL,
                      comment="Used to determine if the Connectome studies "
                              "will be added to the pipeline or not.",
                      values=["Off", "On"],
                      wkf_switch = True)

        self.page.add(label="Phenotype file", 
                     control=control.COMBO_BOX,
                     name="connectome.phenotype.file", 
                     type=dtype.STR, 
                     values="",
                     comment="Path to a CSV file containing the phenotypic "
                             "information.")
        
        self.page.add(label="Phenotype Participant Column",
                      control=control.TEXT_BOX,
                      name="connectome.phenotype.participant_column",
                      type=dtype.STR,
                      comment="Name of the participants column in your "
                              "regressor file.",
                      values="")

        self.page.add(label="Phenotype Variable of Interest Column", 
                     control=control.TEXT_BOX, 
                     name="connectome.phenotype.variable_column", 
                     type=dtype.STR,
                     values="",
                     comment="Column from the CSV file indicating factor "
                             "variables. Only binary variables are valid.")

        self.page.add(label="Phenotype Stratification Column", 
                     control=control.TEXT_BOX, 
                     name="connectome.phenotype.stratification_column", 
                     type=dtype.STR,
                     values="",
                     comment="Column from the CSV file indicating factor "
                             "variables. Only binary variables are valid.")

        self.page.add(label="Folds",
                      control=control.INT_CTRL,
                      name="connectome.folds",
                      type=dtype.NUM,
                      comment="Number of folds to execute cross-validation.",
                      values=5)

        self.page.add(label="Regions of Interest",
                      control=control.LISTBOX_COMBO,
                      name='connectome.roi',
                      comment="",
                      size=(400, 100),
                      combo_type=RoiEditor,
                      type=dtype.LOBJ,
                      type_metadata={
                          "class": RoiOption
                      })

        self.page.add(label="Connectivity Measures",
                      control=control.LISTBOX_COMBO,
                      name='connectome.connectivity',
                      comment="",
                      size=(400, 100),
                      combo_type=ConnectivityEditor,
                      type=dtype.LOBJ,
                      type_metadata={
                          "class": ConnectivityOption
                      })

        self.page.add(label="Classifiers",
                      control=control.LISTBOX_COMBO,
                      name='connectome.classifier',
                      comment="",
                      size=(400, 100),
                      combo_type=ClassifierEditor,
                      type=dtype.LOBJ,
                      type_metadata={
                          "class": ClassifierOption
                      })

        self.page.set_sizer()
        parent.get_page_list().append(self)
        
    def get_counter(self):
        return self.counter
