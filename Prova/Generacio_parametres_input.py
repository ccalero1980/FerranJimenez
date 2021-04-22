import numpy as np
import sys
# PARAMETRES: Introdueixo la frequencia i la magnitud del camp magnetic exteriorment

freq_field = float(sys.argv[1]) #Hz
magB = float(sys.argv[2]) # Tesla
B = magB
Bz = magB # Tesla
# ----------------------------------------

long_filament = 10.e-6 # m
radi_filament = 2.e-7  #m
densitat = 8000. # kg/m^3
massa_filament = densitat*3.1416*radi_filament**2.*long_filament # kg
viscositat = 0.001 #Pa*s
gamma_filament = 2.*3.1416*viscositat*long_filament/(np.log(long_filament/radi_filament) + np.log(2)-0.5)
bend_stiffness = 5.e-20 # J*m (ordre de magnitud, en el cas estudiat pel Gauger es de 4.5e-22)

mu0 = 4*3.1416*1.e-7 # T^2*m^3/J = T^2*m*s^2/kg

Ncycles = 10 # Nombre de cicles del camp magnetic a simular
moment_magnetic_filament = 3.72e-14 # J/T
gravetat = 10. #m/s^2

t0 = 1.e-6 # segons
tinercial = massa_filament/gamma_filament

oseen = 0
blake = 1

fpars = open('parameters.dat','w')

fpars.write( '\n\t PARAMETRES DEL PROBLEMA:\n')
fpars.write( '\nlong_filament (m) = %le\nmassa_filament (kg) = %le\nt0 = %le\n'%(long_filament, massa_filament, t0))
fpars.write( 'viscositat (en kg/(m*s)) = %le\ngamma_filament (en kg/s) = %le\n'%(viscositat, gamma_filament))
fpars.write( 'bend_stiffness (en m^3*kg/s^2) = %le\n'%bend_stiffness)
fpars.write( 'Bfield (T) = %le\nfreq_field (1/s) = %le\n'%(B, freq_field))
fpars.write( 'moment_magnetic_filament (en kg*m^2/(Tesla s^2) = %le\n'%(moment_magnetic_filament))
fpars.write( 'tinercial = %le (segons)\n'%tinercial)


print '\nlong_filament (m) = %le, massa_filament (kg) = %le, t0 = %le'%(long_filament, massa_filament, t0)
print 'viscositat (en kg/(m*s)) = %le, gamma_filament (en kg/s) = %le'%(viscositat, gamma_filament)
print 'bend_stiffness (en m^3*kg/s^2) = %le'%bend_stiffness
print 'Bfield (T) = %le, freq_field (1/s) = %le'%(B, freq_field)
print 'moment_magnetic_filament (en kg*m^2/(Tesla s^2) = %le'%(moment_magnetic_filament)
print 'tinercial = %le (segons)'%tinercial

# CANVI UNITATS:
# Unitats del problema: long_filament, massa_filament, t0 = 10^-7s (de manera que viscositat = 1 en unitats del problema).

radi_filament = radi_filament/long_filament  # en long_filament
viscositat = viscositat*(long_filament/massa_filament)*t0 # en massa_filament/(long_filament*t0)
gamma_filament =2.*3.1416*viscositat/(np.log(1./radi_filament) + np.log(2)-0.5)
bend_stiffness = (bend_stiffness/(long_filament**3.*massa_filament))*t0**2 #  en long_filament**3 * massa_filament/t0^2

B = B # Tesla
freq_field = freq_field*t0 #1/t0
moment_magnetic_filament = (moment_magnetic_filament/(massa_filament*long_filament**2.))*t0**2. # en massa_filament*long_filament**2/(t0^2*Tesla)
mu0 = mu0*(massa_filament/(long_filament*t0**2)) # mu0 en T^2*long_filament*t0^2/massa_filament
gravetat = gravetat*t0**2/long_filament

massa_filament = massa_filament/massa_filament
long_filament = long_filament/long_filament

fpars.write( '\n\nCANVI UNITATS:\n')
fpars.write( 'long_filament (lf) = %f\nradi_filament (en lf)= %f\nmassa_filament (mf) = %f\n'%(long_filament, radi_filament, massa_filament))
fpars.write( 'viscositat (en mf/(lf*t0)) = %f\ngamma_filament (en mf/t0) = %f\n'%(viscositat, gamma_filament))
fpars.write( 'bend_stiffness (en lf^3*mf/t0^2) = %f\n'%bend_stiffness)
fpars.write( 'Bfield (T) = %f\nfreq_field (1/t0) = %f\n'%(B, freq_field))
fpars.write( 'moment_magnetic_filament (en mf*lf^2/(Tesla t0^2) = %f\n'%(moment_magnetic_filament))

print '\nCANVI UNITATS:'
print 'long_filament (lf) = %f, radi_filament (en lf)= %f, massa_filament (mf) = %f'%(long_filament, radi_filament, massa_filament)
print 'viscositat (en mf/(lf*t0)) = %f, gamma_filament (en mf/t0) = %f'%(viscositat, gamma_filament)
print 'bend_stiffness (en lf^3*mf/t0^2) = %f'%bend_stiffness
print 'Bfield (T) = %f, freq_field (1/t0) = %f'%(B, freq_field)
print 'moment_magnetic_filament (en mf*lf^2/(Tesla t0^2) = %f'%(moment_magnetic_filament)
# QUANTITATS BASIQUES

tinercial = (massa_filament/gamma_filament) # en t0

print 'tinercial = %.8f (t0)'%tinercial
if tinercial > 0.1*(1./freq_field): print 'ERROR tinercial'


dt = tinercial/100.
Tsim = Ncycles*(2.*np.pi/(freq_field))
Ntimesteps = int(Tsim/dt)

fpars.write( '\n\t SIMULACIO:\n')
fpars.write( 'timestep = %le (t0)\nSimulation_Time = %le (t0)\nNtimesteps = %i\n'%(dt, Tsim, Ntimesteps))
if oseen==1:fpars.write('OSEEN\n')
elif blake==1:fpars.write('BLAKE\n')
print 'timestep = %le (t0), Simulation_Time = %le (t0), Ntimesteps = %i'%(dt, Tsim, Ntimesteps)
#
if freq_field <100*t0:
  steps_in_cycle = 10000
elif freq_field >= 100*t0 and freq_field<=1000*t0:
  steps_in_cycle = 1000
elif freq_field >1000*t0:
  steps_in_cycle = 100
  
number_cycles = int(float(Ntimesteps)/float(steps_in_cycle))


fout = open('init', 'w')
fout.write('./conf_inicial_straight10  : file name\n')
fout.write('./       : output dir\n')
fout.write('%i	     : number of cycles\n'%number_cycles)
#fout.write('%i	     : number of cycles\n'%1)
fout.write('%f       : timestep\n'%dt)
fout.write('%f       : bending modulus\n'%bend_stiffness)
fout.write('%f       : mass filament\n'%massa_filament)
fout.write('%f       : viscosity\n'%viscositat)
fout.write('%f       : gamma filament\n'%gamma_filament)
fout.write('0.3           : hydr radius (in bond lengths)\n')
fout.write('1.0e-8        : constraint tolerance\n')
fout.write('%i	     : steps per cycle\n'%steps_in_cycle)
fout.write('%f       : Bx\n'%B)
fout.write('%f       : Bz\n'%Bz)
fout.write('%f       : freq_field \n'%freq_field)
fout.write('%f       : moment_magnetic_filament\n'%moment_magnetic_filament)
fout.write('%f       : mu0 permeability vacuum \n'%mu0)
fout.write('%f       : gravetat \n'%gravetat)
fout.write('%i       : Oseen\n'%oseen)
fout.write('%i       : Blake\n'%blake)








fout.close()
