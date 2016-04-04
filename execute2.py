
from technologies import silicon_photonics
from ipkiss3 import all as i3
from picazzo3.traces.wire_wg import WireWaveguideTemplate
from spiralmzi import SpiralMZI
import numpy as np
import pylab as plt

radius = 20.0

start = 1.54
stop = 1.56
samples = 500
wavelengths = np.linspace(start, stop, samples)
loss_dB_cm = 3
wg = WireWaveguideTemplate()
wg_layout=wg.Layout(core_width=0.45)
my_mzi=SpiralMZI(loss_dB_cm=loss_dB_cm,radius=radius,dc_coupler_length=20.6, dc_gap=0.2, wg_template=wg, wg_dc_template=wg)
my_mzi_layout = my_mzi.Layout()
my_mzi_layout.visualize()

reme_model = my_mzi.RemeSimulation(start=start,stop=stop)
reme_model.calculation()

n_eff_straight = my_mzi.RemeSimulation.view.fstraight
plt.plot(wavelengths, [n_eff_straight(w).real for w in wavelengths], label='Straight')
plt.xlim(start, stop)
plt.xlabel("Wavelength ($\mu m$)")
plt.ylabel("Neff")
plt.legend()
plt.show()

reme_model.straight_guide.plot_mode(0, 'Ey', 'Real')
reme_model.coupler_dc.plot_transfer_matrix()

my_mzi_cm = my_mzi.CapheModel()
wavelengths = np.linspace(start,stop,samples)
s_matrix = my_mzi_cm.get_smatrix(wavelengths=wavelengths)
#engine = i3.CapheFrequencyEngine()
#simulation = engine.SMatrixSimulation(model=my_mzi_cm,
#                                     wavelengths=wavelengths)
#simulation.run()
#s_matrix = simulation.monitors['s_matrix']
#print s_matrix
plt.plot(wavelengths, np.abs(s_matrix['in_1', 'out_1'])**2, label='Straight transmission')
plt.plot(wavelengths, np.abs(s_matrix['in_1', 'out_2'])**2, label='Cross transmission')
plt.title("my mzi")
plt.xlabel("wavelength ($\mu m$)")
plt.ylabel("Power transmission")
plt.legend()
plt.show()
