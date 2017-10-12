#5.0/surfrace
cd fillin
from LOADRUNREQs import * 
cd /scratch/midway/jriback/runupside/PNt_ALLH_0_MP0p7_0_sim.h5tmp
for x in os.listdir('/fillin/runupside/PNt_ALLH_0_MP0p7_0_sim.h5tmpFOXSall'):
    if (x[-3:]=='int')and(x[:2]=='xx'):
        sp.check_output('printf \'1\n%s\n1.4\n1\n\' | /fillin/surfrace/5.0/surfrace &&' %(x[:-6]+'.pdb'))
        print x
        break
