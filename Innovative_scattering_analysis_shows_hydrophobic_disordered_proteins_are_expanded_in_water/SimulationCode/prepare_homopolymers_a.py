from LOADRUNREQs import *

proti=1
#MAke Average PNt RAMAs
#jobs = dict()
AVERAGERAMA={}
proteins=[('PNt','AverageRama'),('polyA','AverageRama'),('polyG','AverageRama')]
numx=33*len(proteins)
print (numx,33)
figure(figsize=(33,numx))
textppii='PPII'
textbeta='BETA'
textalphar='APH(R)'
textalphal='APH(L)'
textgamma='GAMMA'
for nm,TYPE in proteins:
    #prefix='%.2f_%.2f_%.2f' % (opts[0],opts[1],opts[2])
        GROUPTORUN=[]
        run=True
        for PPII in [0]:#range(0,11):
            for LJDFIRE in [0]:#range(0,11):
                prefix='%s_%s' % (nm,TYPE)

                config=run_dir+'%s_%s_%s_%s_sim.h5'%(nm,TYPE,PPII,LJDFIRE)

                #ramafile=run_dir+'%s_%s.rama.pkl'%(nmi,TYPE)

                if os.path.exists(config):
                    print config

                    #Load TrueRama

                    file2=tb.open_file(config,'r')
                    AVERAGERAMA[nm]=0.
                    size=file2.root.input.potential.rama_map_pot.rama_pot.shape[0]
                    print nm,size
                    cnt=0.0
                    for x in range(1,size-1):
                        AVERAGERAMA[nm]+=(np.exp(-file2.root.input.potential.rama_map_pot.rama_pot[x])/np.exp(-file2.root.input.potential.rama_map_pot.rama_pot[x]).sum())
                        #print np.exp(-file2.root.input.potential.rama_map_pot.rama_pot[x]).sum()
                        cnt+=1.0
                    print size,cnt,AVERAGERAMA[nm].sum()
                    file2.close()
                    
                    tmpA=AVERAGERAMA[nm]/cnt
                    tmpA=tmpA/(tmpA.sum())
                    
                    if True:
                        ppii=0.0
                        for i in range(72):
                            if not ((i>45)|(i<=15)):continue
                            for j in range(72):
                                if not (15<j<=35):continue
                                ppii+=tmpA[j,i]
                        beta=0.0
                        for i in range(72):
                            if not ((i>45)|(i<=15)):continue
                            for j in range(72):
                                if not (15>=j):continue
                                beta+=tmpA[j,i]
                        alpha=0.0
                        for i in range(72):
                            if not (15<i<=45):continue
                            for j in range(72):
                                if not (j<=35):continue
                                alpha+=tmpA[j,i]
                        alphal=0.0
                        for i in range(72):
                            if not (25<i<=55):continue
                            for j in range(72):
                                if not (j>35):continue
                                alphal+=tmpA[j,i]
                        try:
                            textppii+='\t%.2f' % (ppii)
                            textbeta+='\t%.2f' % (beta)
                            textalphar+='\t%.2f' % (alpha)
                            textalphal+='\t%.2f' % (alphal)
                            textgamma+='\t%.2f' % (1-alpha-ppii-beta-alphal)
                        except tb.NoSuchNodeError:
                            print 'FAIL'
                    
                    
                    AVERAGERAMA[nm]=-np.log(AVERAGERAMA[nm]/cnt)
                    tmpRAMA=AVERAGERAMA[nm]
                    if False:
                        ppii=0.0
                        for i in range(72):
                            if not ((i>45)|(i<=15)):continue
                            for j in range(72):
                                if not (15<j<=35):continue
                                tmpRAMA[j,i]+=15
                    if False:
                        beta=0.0
                        for i in range(72):
                            if not ((i>45)|(i<=15)):continue
                            for j in range(72):
                                if not (15>=j):continue
                                tmpRAMA[j,i]+=15
                    if False:
                        alpha=0.0
                        for i in range(72):
                            if not (15<i<=45):continue
                            for j in range(72):
                                if not (j<=35):continue
                                tmpRAMA[j,i]+=15
                    if False:
                        for i in range(72):
                            if not (25<i<=55):continue
                            for j in range(72):
                                if not (j>35):continue
                                tmpRAMA[j,i]=-15
                    
                    subplot(3,3,proti)
                    proti+=1
                    d2=np.swapaxes(-tmpRAMA,0,1)#np.exp(-t.root.input.potential.rama_map_pot.rama_pot[5]),0,1)
                    
                    imshow(-d2+d2.max(), origin='lower left', extent=(-180,180,-180,180),vmin=0.0,vmax=5.0)
                    
                    continue
                    #run=False
                    continue


                fastafile = run_dir+'%s.fasta'% nm
                
                print 'STARTING FINAL',nm,TYPE,PPII,LJDFIRE

                sidechaini='PNt_%s.h5' % LJDFIRE
                if LJDFIRE==None:
                    sidechaini=None
                
                upside_config('%s.fasta' % nm, config, 1,
                      hbond=False, reference_rama=None, 
                      sidechain=sidechaini,rama_pot='rama_libraries.h5',rama_pot_product=True,#sidechain_scale=MJi, inverse_scale=opts[1], inverse_radius_scale=opts[2],
                      rotamer=False)
print textppii
print textbeta
print textalphar
print textalphal
print textgamma
