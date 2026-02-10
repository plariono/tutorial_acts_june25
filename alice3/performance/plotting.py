import acts
import acts.examples

#Residual Plots

resPlotToolConfig = acts.examples.root.ResPlotToolConfig()
binning = resPlotToolConfig.varBinning
binning["Eta"] = acts.examples.root.AxisVariant.regular(80, -4,    4,     "#eta")
binning["Pt"] =  acts.examples.root.AxisVariant.regular(1000, -0,    50,     "pT [GeV/c]")
binning["Residual_phi"] =  acts.examples.root.AxisVariant.regular(200, -0.02,    0.02,     "r_{#phi} [rad]")
binning["Residual_theta"] =  acts.examples.root.AxisVariant.regular(200, -0.02,    0.02,     "r_{#theta} [rad]")
binning["Residual_qop"] =  acts.examples.root.AxisVariant.regular(200, -0.2,    0.2,     "r_{q/p} [c/GeV]")
resPlotToolConfig.varBinning = binning

# Duplication Plots
duplicationPlotToolConfig = acts.examples.root.DuplicationPlotToolConfig()
duplicationPlotToolConfig.varBinning = {
    "Eta": acts.examples.root.AxisVariant.regular(80, -4,    4,     "#eta"),
    "Phi": acts.examples.root.AxisVariant.regular(100, -3.15, 3.15, "#phi"),
    "Pt":  acts.examples.root.AxisVariant.regular(1000,  0,    50,   "pT [GeV/c]"),
"Num": acts.examples.root.AxisVariant.regular(30, -0.5,  29.5,  "N"),
}

# Efficiency Plots
effPlotToolConfig = acts.examples.root.EffPlotToolConfig()
binning = effPlotToolConfig.varBinning
binning["Eta"] = acts.examples.root.AxisVariant.regular(80, -4,    4,     "#eta")
binning["Pt"]  = acts.examples.root.AxisVariant.regular(1000,  0,    50,   "pT [GeV/c]")
effPlotToolConfig.varBinning = binning
effPlotToolConfig.minTruthPt = 0.15

# Fake rate Plots
fakePlotToolConfig = acts.examples.root.FakePlotToolConfig()
binning = fakePlotToolConfig.varBinning
binning["Eta"] = acts.examples.root.AxisVariant.regular(80, -4,    4,     "#eta")
binning["Pt"] = acts.examples.root.AxisVariant.regular(1000,  0,    50,   "pT [GeV/c]")
fakePlotToolConfig.varBinning = binning

# TrackQuality Plots
trackQualityPlotToolConfig = acts.examples.root.TrackQualityPlotToolConfig()
binning = trackQualityPlotToolConfig.varBinning
binning["Eta"] = acts.examples.root.AxisVariant.regular(80, -4,    4,     "#eta")
binning["Pt"] = acts.examples.root.AxisVariant.regular(1000,  0,    50,   "pT [GeV/c]")
trackQualityPlotToolConfig.varBinning = binning


# Track Summary Plots
trackSummaryPlotToolConfig = acts.examples.root.TrackSummaryPlotToolConfig()
binning = trackSummaryPlotToolConfig.varBinning
binning["Eta"] = acts.examples.root.AxisVariant.regular(80, -4,    4,     "#eta")
binning["Pt"] = acts.examples.root.AxisVariant.regular(1000,  0,    50,   "pT [GeV/c]")
trackSummaryPlotToolConfig.varBinning = binning
