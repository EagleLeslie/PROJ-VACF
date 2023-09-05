import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.fft import fft, fftfreq

time_001, real_001, imag_001 = np.loadtxt(fname = './fort_001.22',unpack=True)
time_010, real_010, imag_010 = np.loadtxt(fname = './fort_010.22',unpack=True)
time_100, real_100, imag_100 = np.loadtxt(fname = './fort_100.22',unpack=True)
tt, v1, v2, v3 = np.loadtxt(fname = './vacfout',unpack=True)

wavenumber = 1E+15/2.998E10
timestep = 1

hhfft_100 = scipy.fft.rfft(real_100[0:1000])
wfft_100 = scipy.fft.rfftfreq(np.size(real_100[0:1000]),d=timestep) * wavenumber

hhfft_010 = scipy.fft.rfft(real_010[0:1000])
wfft_010 = scipy.fft.rfftfreq(np.size(real_010[0:1000]),d=timestep) * wavenumber

hhfft_001 = scipy.fft.rfft(real_001[0:1000])
wfft_001 = scipy.fft.rfftfreq(np.size(real_001[0:1000]),d=timestep) * wavenumber

# hhfft_vacf = scipy.fft.rfft(v1[0:500])
# wfft_vacf = scipy.fft.rfftfreq(np.size(v1[0:500]),d=timestep) * wavenumber

# plt.plot(tt,v1,label='Fe')
plt.plot(tt,v2,color='darkorange',label='H')
# plt.plot(tt,v3,label='O')
P
plt.plot(time_001,real_001,"--",color='red',label='q=(0,0,1)')
plt.plot(time_010,real_010,"--",color='forestgreen',label='q=(0,1,0)')
plt.plot(time_100,real_100,"--",color='purple',label='q=(1,0,0)')
plt.axhline(y=0,color='black',zorder=0)
plt.xlabel('Timestep (fs)')
plt.ylabel('Velocity Autocorrleation Function')
plt.xlim(0,300)
plt.legend()
plt.show()

# plt.plot(time2,imag2)
# plt.show()

plt.plot(wfft_100, hhfft_100,label='q=(1,0,0)')
plt.plot(wfft_010, hhfft_010,label='q=(0,1,0)')
plt.plot(wfft_001, hhfft_001,label='q=(0,0,1)')
plt.xlabel('Wavenumber (cm$^{-1}$)')
plt.ylabel('Intensity (arbitrary units)')
plt.xlim(0,3000)
plt.legend()
plt.show()