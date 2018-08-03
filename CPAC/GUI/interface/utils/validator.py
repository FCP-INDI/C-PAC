import wx
import string

class NotEmptyValidator(wx.PyValidator): 

    def __init__(self):        
        wx.PyValidator.__init__(self)
        
    def Clone(self):
        """
        Note that every validator must implement the Clone() method.
        """
        return NotEmptyValidator()

    def Validate(self, win):  
        
        textCtrl = self.GetWindow()
        
        text = textCtrl.GetValue()
        if len(text) == 0:
            wx.MessageBox("This field must contain some text!", "Error")
            textCtrl.SetBackgroundColour("pink")
            textCtrl.SetFocus()
            textCtrl.Refresh()
            return False
        else:
            textCtrl.SetBackgroundColour(
                wx.SystemSettings_GetColour(wx.SYS_COLOUR_WINDOW))
            textCtrl.Refresh()
            return True

    def TransferToWindow(self):
            return True
    
    def TransferFromWindow(self):
            return True
            
        
class CharValidator(wx.PyValidator):
    def __init__(self, flag):
        wx.PyValidator.__init__(self)
        self.flag = flag
        self.Bind(wx.EVT_CHAR, self.OnChar)

    def Clone(self):
        return CharValidator(self.flag)
    
    def Validate(self, win):
        return True
    
    def TransferToWindow(self):
        return True
    
    def TransferFromWindow(self):
        return True
    
    def OnChar(self, evt):
        key = chr(evt.GetKeyCode())
        if self.flag == "no-alpha" and key in string.letters:
            return
        if self.flag == "no-digit" and key in string.digits:
            return
        evt.Skip()
