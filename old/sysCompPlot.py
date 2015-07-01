from ROOT import *
from array import array
import os

x                = array( 'd', [1.05,   11.75,  15.5,   17.0,   19.0,   17.5    ])

fit_bias_afb     = TGraph(6,x,array( 'd', [0.011,  0.0012, 0.0059, 0.0097, 0.0040, 0.0005  ])  )
resolution_afb   = TGraph(6,x,array( 'd', [0.0011, 0.0016, 0.0006, 0.0014, 0.0014, 0.0013  ])  )
#non_flat_eff_afb = array( 'd', [1.05,11.75,15.5,17.0,19.0,17.5])
bkg_shape_afb    = TGraph(6,x,array( 'd', [0.003,  0.045,  0.010,  0.026,  0.011,  0.007   ])  )
mc_stat_afb      = TGraph(6,x,array( 'd', [0.0015, 0.0012, 0.0009, 0.0007, 0.0003, 0.0039  ])  )
mismodel_afb     = TGraph(6,x,array( 'd', [0.002, 0.0069,  0.0018, 0.0012, 0.0030, 0.0002  ])  )


fit_bias_afbB     = TGraph(6,x,array( 'd', [0.14,  0.0203, 0.0077, 0.0016, 0.0026, 0.0010  ])  )
resolution_afbB   = TGraph(6,x,array( 'd', [0.007, 0.013,  0.010,  0.010,  0.010,  0.0011  ])  )
#non_flat_eff_afbB = array( 'd', [1.05,11.75,15.5,17.0,19.0,17.5])
bkg_shape_afbB    = TGraph(6,x,array( 'd', [0.053,  0.035,  0.026,  0.022,  0.025,  0.017   ])  )
mc_stat_afbB      = TGraph(6,x,array( 'd', [0.0021, 0.0020, 0.0010, 0.0042, 0.0016, 0.0014  ])  )
mismodel_afbB     = TGraph(6,x,array( 'd', [0.0093, 0.0069, 0.0109, 0.0159, 0.0148, 0.0138  ])  )


fit_bias_fL     = TGraph(6,x,array( 'd', [0.0400,  0.0360, 0.0031, 0.0012, 0.0031, 0.0042  ])  )
resolution_fL   = TGraph(6,x,array( 'd', [0.0135,  0.0135, 0.0031, 0.0008, 0.0013, 0.0008  ])  )
#non_flat_eff_fL = array( 'd', [1.05,11.75,15.5,17.0,19.0,17.5])
bkg_shape_fL    = TGraph(6,x,array( 'd', [0.049,   0.034,  0.038,  0.036,  0.031,  0.014   ])  )
mc_stat_fL      = TGraph(6,x,array( 'd', [0.0017, 0.0015,  0.0002, 0.0025, 0.0037, 0.0009  ])  )
mismodel_fL     = TGraph(6,x,array( 'd', [0.044,  0.0027,  0.0046, 0.0043, 0.0017, 0.0046  ])  )


names =   [ "fit bias"  , "resolution",   "bkg shape",   "MC stat",   "Mismodel" ]
afb_grs = [ fit_bias_afb, resolution_afb, bkg_shape_afb, mc_stat_afb, mismodel_afb ]
afbB_grs = [ fit_bias_afbB, resolution_afbB, bkg_shape_afbB, mc_stat_afbB, mismodel_afbB ]
fL_grs = [ fit_bias_fL, resolution_fL, bkg_shape_fL, mc_stat_fL, mismodel_fL ]

varnames = ["afb", "afbB", "fL"]
sys_dict = dict(zip(varnames,[afb_grs, afbB_grs, fL_grs]))

for v in varnames :
	
	c = TCanvas()
	multigr = TMultiGraph(v,"Systematics on "+v)
	leg = TLegend(0.7,0.7,0.99,0.99)

	for i,g in enumerate(sys_dict[v]) :
		g.SetMarkerStyle(22+i)
		g.SetMarkerSize(1)
		g.SetMarkerColor(2+i)
		multigr.Add(g)
		leg.AddEntry(g,names[i],"P")

	multigr.Draw("AP")
	multigr.GetXaxis().SetTitle("#it{q}^{2} [Gev^{2}/#it{c}^{4}]")
	multigr.GetYaxis().SetTitle("Abs. Sys.")
	leg.Draw()
	c.Print(v+"_sys.pdf")
	os.system("sleep 10")





