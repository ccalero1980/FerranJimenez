import sys
numin = int(sys.argv[1])
print numin
nfil = 1
nparts1 = 10

nmons = nparts1

nameout = 'trajectory.xyz'
fout = open(nameout,'w')
fout.write(str(nmons) + '\n\n')
fin = open('conf_inicial_straight10', 'r')
for k in range(nparts1):
  line = fin.readline()
  x=float(line.split()[0])
  y=float(line.split()[1])
  z=float(line.split()[2])
  fout.write('A  %14.8f  %14.8f  %14.8f\n'%(x,y,z))

for n in range(1,numin):
  fout.write(str(nmons) + '\n\n')
  namein = 'out_%i.00'%n
  fin = open(namein, 'r')
  for k in range(nparts1):
    line = fin.readline()
    x=float(line.split()[0])
    y=float(line.split()[1])
    z=float(line.split()[2])
    fout.write('A  %14.8f  %14.8f  %14.8f\n'%(x,y,z))
 
