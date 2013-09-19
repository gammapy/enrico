import ROOT


def RootStyle():
    """Make ROOT plots look less ugly"""
    ROOT.gROOT.SetStyle("Plain")
    ROOT.gStyle.SetPalette(1)
    # Canvas
    ROOT.gStyle.SetCanvasColor(10)
    # Frame
    ROOT.gStyle.SetFrameBorderMode(0)
    ROOT.gStyle.SetFrameFillColor(0)
    # Pad
    ROOT.gStyle.SetPadBorderMode(0)
    ROOT.gStyle.SetPadColor(0)
    ROOT.gStyle.SetPadTopMargin(0.07)
    ROOT.gStyle.SetPadLeftMargin(0.13)
    ROOT.gStyle.SetPadRightMargin(0.11)
    ROOT.gStyle.SetPadBottomMargin(0.1)
    ROOT.gStyle.SetPadTickX(1)
    # Make ticks be on all 4 sides
    ROOT.gStyle.SetPadTickY(1)
    # histogram
    ROOT.gStyle.SetHistFillStyle(0)
    ROOT.gStyle.SetOptTitle(0)
    # histogram title
    ROOT.gStyle.SetTitleSize(0.22)
    ROOT.gStyle.SetTitleFontSize(2)
    ROOT.gStyle.SetTitleFont(42)
    ROOT.gStyle.SetTitleFont(62, "xyz")
    ROOT.gStyle.SetTitleYOffset(1.0)
    ROOT.gStyle.SetTitleXOffset(1.0)
    ROOT.gStyle.SetTitleXSize(0.04)
    ROOT.gStyle.SetTitleYSize(0.04)
    ROOT.gStyle.SetTitleX(.15)
    ROOT.gStyle.SetTitleY(.98)
    ROOT.gStyle.SetTitleW(.70)
    ROOT.gStyle.SetTitleH(.05)
    # statistics box
    ROOT.gStyle.SetStatFont(42)
    ROOT.gStyle.SetStatX(.91)
    ROOT.gStyle.SetStatY(.90)
    ROOT.gStyle.SetStatW(.15)
    ROOT.gStyle.SetStatH(.15)
    # axis labels
    ROOT.gStyle.SetLabelFont(42, "xyz")
    ROOT.gStyle.SetLabelSize(0.035, "xyz")
    ROOT.gStyle.SetGridColor(16)
    # legend
    ROOT.gStyle.SetLegendBorderSize(0)
