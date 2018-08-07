import wx
import wx.html
from ..utils.generic_class import GenericClass
from ..utils.constants import control, dtype
from ..utils.custom_control import ListBoxComboEditor
import os
import pkg_resources as p


class RoiOption:
    
    labels = {
        'dict_learning': 'Dictionary Learning',
        'k_means': 'k-Means',
        'ward': 'Ward',
        'group_ica': 'Group ICA',
    }

    def __init__(self, type, parameters):
        assert type in self.labels
        self.type = type
        self.parameters = parameters

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return '%s [%s]' % (
            self.labels[self.type],
            ";".join(['%s: %s' % (k, v) for k, v in self.parameters.items()])
        )


class RoiEditor(ListBoxComboEditor):

    def __init__(self, parent):
        super(RoiEditor, self).__init__(
            parent=parent, title="Region of Interest Selector",
            size=(300, 300)
        )

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

    def __iter__(self):
        return iter([])

    
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
                      values="",
                      comment="",
                      size=(400, 100),
                      combo_type=RoiEditor,
                      type=dtype.LOBJ,
                      type_metadata={
                          "class": RoiOption
                      })

        self.page.set_sizer()
        parent.get_page_list().append(self)
        
    def get_counter(self):
        return self.counter
