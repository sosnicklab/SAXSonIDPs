import os
import numpy as np
tmpstrng=''
for interact in range(30):
    ASAz=[]
    for rep in range(5):
        print interact,rep
        ASAs=[]
        cnt=0
        tmpstrng+='%s\t%s' % (interact,rep)
        for FILE in list((x for x in os.listdir('/fillin/runupside/PNt_ALLH_%s_MP0p7_%s_sim.h5tmpFOXSall' % (interact,rep)) if (x[-3:]=='txt')and(x[:2]=='xx')and(len(x.split('_'))==2))):
            cnt+=1
            with open('/fillin/runupside/PNt_ALLH_%s_MP0p7_%s_sim.h5tmpFOXSall/%s' % (interact,rep,FILE),'r') as F:
                ASA=0.0
                for LINE in F.readlines():
                    ASA+=float(LINE.strip().split(' ')[-1])
                if ASA<1000.:continue
                if ASA>45000.:continue
                ASAs.append(ASA)
                tmpstrng+='\t%.2f' % ASA
            #print FILE,ASA
        tmpstrng+='\n'
        print "HERE=",interact,rep,np.mean(ASAs),cnt
fileasas=open('/fillin/runupside/PNt_ALLH_MP0p7_FOXSall_ASA.txt','w')
fileasas.write(tmpstrng)
fileasas.close()
