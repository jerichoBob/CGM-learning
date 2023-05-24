 # no_minimize.py
# you gotta run this with pythonw...
#   which means you've' done this:
#     conda install -c conda-forge python.app

import wx

app = wx.App()
frame = wx.Frame(None, style=wx.MAXIMIZE_BOX | wx.RESIZE_BORDER
	| wx.SYSTEM_MENU | wx.CAPTION |	 wx.CLOSE_BOX)
frame.Show(True)

app.MainLoop()