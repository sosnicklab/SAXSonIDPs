###THIS IS FOR FOXS!!!!!
#for x in {0..29}; do echo $x; less /fillin/PNt_ALLH_"$x"_MP0p7_4_sim.h5_vtfetcFOXSall.output |grep foxs|wc -l; done
from LOADRUNREQs import *

parser = argparse.ArgumentParser(description='Prepare input file',
            usage='use "%(prog)s --help" for more information')
parser.add_argument('--configA', required=True,
            help='config file name before param number')
parser.add_argument('--configB', required=True,
            help='config file name after param number before rep number _sim.h5')
parser.add_argument('--Params', required=True,
            help='config file name without rep.h5')
parser.add_argument('--configC', required=False,default='',
            help='special options')
args = parser.parse_args()

PERFECT=True

if args.configC=="":
    print 'WRONG PROCESS PYTHON FORMULA'
    exit()
if args.configC=="ALA":
    print 'WRONG PROCESS PYTHON FORMULA'
    exit()

for MPDFIRE in range(0,int(args.Params)):
    for rep in range(0,5):
        config='%s_%s_%s_%s_sim.h5' % (args.configA,MPDFIRE,args.configB,rep)
        print config
        os.popen('rm %stmp%s/xx50000*' % (config,args.configC))
        if len(list(( x for x in os.listdir('%stmp%s' % (config,args.configC)) if (x[-3:]=='dat')and(x[:2]=='xx'))))!=500:
            print 'only found %s out of 500' % len(list(( x for x in os.listdir('%stmp%s' % (config,args.configC)) if (x[-3:]=='dat')and(x[:2]=='xx'))))
            PERFECT=False
            #print 'FAILED EXITING'
            #sys.exit()
        print '\tFoundAll'
fileScattering=''
if PERFECT:
    print 'All found---Moving to produce data\n\n'
    print 'saving to file %s_%s_SCATTERING%s' % (args.configA,args.configB,args.configC)
    fileScattering=open('%s_%s_SCATTERING%s' % (args.configA,args.configB,args.configC),'w')
else :
    print 'MISSING SOME---Moving to produce data\n\n'
    print 'saving to file %s_%s_SCATTERING_PARTIAL%s' % (args.configA,args.configB,args.configC)
    fileScattering=open('%s_%s_SCATTERING_PARTIAL%s' % (args.configA,args.configB,args.configC),'w')
Iqi='{'#'----------\n{'
for numEng in range(0,int(args.Params)):
    print numEng
    Iqi+='{'
    configi='%s_%s_%s_%s_sim.h5' % (args.configA,numEng,args.configB,'%s')
    first=True
    for rep in range(0,5):
        config=configi % (rep)
        if not os.path.exists('%stmp%s' % (config,args.configC)):continue
        if len(list(x for x in os.listdir('%stmp%s' % (config,args.configC)) if (x[-3:]=='dat')and(x[:2]=='xx')))==0:continue
        qs=[]
        Is=[]
        #Rgs=np.array(list(float(x.split(':')[2][1:-1]) for x in (os.popen('grep \'Radius of Gyration\' %stmp%s/xx*.log' % (config,args.configC)).readlines())))
        with tb.open_file(config) as t:
            if ((t.root.output.pos[50:,0,1::3]).sum(axis=-1).mean(axis=-1)**2).mean()>0.001:
                print 'PROBLEM WITH CENTERING!'
            Rgs=np.sqrt((t.root.output.pos[50:,0,1::3]**2).sum(axis=-1).mean())
        cnt=0.
        for x in list(x for x in os.listdir('%stmp%s' % (config,args.configC)) if (x[-3:]=='dat')and(x[:2]=='xx')):
            with open('%stmp%s/%s' % (config,args.configC,x)) as filei:
                filedaz=filei.readlines()[2:]
                tmpI=[]
                for x in filedaz:
                    dazi=x.split(' ')
                    #if cnt==0:qs.append(parseE(dazi[0]))
                    tmpI.append(float(dazi[4]))
                Is.append(tmpI)
            cnt+=1.
        #print qs
        norm=np.array(Is).mean(axis=0)[0]
        if not first:
            Iqi+=',\n'
        else:
            first=False
        Iqi+='{\"%s\",%.3f}' % (str(np.array(Is).mean(axis=0)/norm*1000)[1:-1],Rgs)#np.sqrt((Rgs**2).mean()))
    if numEng!=(int(args.Params)-1):
        Iqi+='},\n'
Iqi+='}'
fileScattering.write(Iqi)
fileScattering.close()

print 'saving to file %s_%s_Rij' % (args.configA,args.configB)
fileRij=open('%s_%s_Rij' % (args.configA,args.configB),'w')
Riji=''#'PNt_HP\tRij'
for numEng in range(0,int(args.Params)):
    print numEng
    config='%s_%s_%s_%s_sim.h5' % (args.configA,numEng,args.configB,'%s')
    Riji+='\n%s\n' % numEng
    POSs=[]
    SIZE=''
    for rep in range(0,5):
            try:
                with tb.open_file(config % (rep)) as t:
                    POSs.append(t.root.output.pos[50:,0,1::3])
                if rep==0:
                    SIZE=len(POSs[0][0])
            except tb.NoSuchNodeError:
                print 'NOT STARTED YET'
            except IOError:
                print 'NOT STARTED YET'
    SIZE=len(POSs[0][0])
    for d in [2, 3, 4, 5, 6, 8, 10, 12, 15, 18, 22, 27, 33, 40, 49, 59, 72, 87, 106, 128, 155, 188, 227, 275, 333, 403, 487, 589, 713]:
        if d>SIZE:continue
        meanz=[]
        #stdz=[]
        for POS in POSs:
            meanz.append(np.sqrt(np.array(list(((POS[:,(x+d)]-POS[:,x])**2).sum(axis=-1).mean() for x in range(0,SIZE-d))).mean()))
            #stdz.append(np.sqrt(np.array(list(np.sqrt(((POS[:,(x+d)]-POS[:,x])**2).sum(axis=-1)).std()**2 for x in range(0,SIZE-d))).sum()/500.))
        Riji+='%s %.5f %.5f\n' % (d,np.array(meanz).mean(),np.array(meanz).std()/np.sqrt(len(meanz)))
    Riji+='\n'
fileRij.write(Riji)
fileRij.close()



#Code to analyze everything and save it to document... also delete all (crysol) files after.
