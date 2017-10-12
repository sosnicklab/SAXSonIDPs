from LOADRUNREQs import *

proti=1
#MAke Average PNt RAMAs
#jobs = dict()
AVERAGERAMA={}
proteins=[('polyG','AverageRama'),('polyA','AverageRama'),('PNt','AverageRama')]
RAMATYPE={}
for nm,TYPE in proteins:
        RAMATYPE[nm]=TYPE
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
                    
                    AVERAGERAMA[nm]=-np.log(AVERAGERAMA[nm]/cnt)
                    
                    continue


                fastafile = run_dir+'%s.fasta'% nm
                
                print 'ERROR DOES NOT EXIST RUN fileA',nm,TYPE,PPII,LJDFIRE
                exit()
                
proteins=[('A','polyG'),('A','polyA'),('A','PNt')]
cnt=1
numx=5*11*len(proteins)
print (numx,33)
figure(figsize=(33,numx))
for nm,homorama in proteins:
    clear_output()
    for rep in [0]:#range(0,5):
        GROUPTORUN=[]
        prefix='%s_%s' % (homorama,rep)
        if RAMATYPE[homorama]=='AverageRamaB':prefix+='b'
        #check ramafile does not exist
        ramafile=run_dir+'%s_%s.h5.rama.pkl'%(nm,prefix)
        if os.path.exists(ramafile):
            print ramafile,'EXITSTS'
            ramax = run_dir+'%s_%s.h5.rama.pkl'%(nm,prefix)
            ALLRAMAi=np.swapaxes(cp.load(open(ramax))[3],0,1)
            AVERAGERAMAfix=AVERAGERAMA[homorama]+np.log(ALLRAMAi/sum(sum(ALLRAMAi))+0.000000001)
            d2A=np.swapaxes(-AVERAGERAMA[homorama],0,1)
            d2B=np.swapaxes(-np.log(ALLRAMAi/sum(sum(ALLRAMAi))+0.000000001),0,1)
            d2C=np.swapaxes(-AVERAGERAMAfix,0,1)#np.exp(-t.root.input.potential.rama_map_pot.rama_pot[5]),0,1)
            subplot(5*len(proteins),3,cnt)
            cnt+=1
            imshow(-d2A+d2A.max(), origin='lower left', extent=(-180,180,-180,180),vmin=0.0,vmax=5.0)
            subplot(5*len(proteins),3,cnt)
            cnt+=1
            imshow(-d2B+d2B.max(), origin='lower left', extent=(-180,180,-180,180),vmin=0.0,vmax=5.0)
            subplot(5*len(proteins),3,cnt)
            cnt+=1
            imshow(-d2C+d2C.max(), origin='lower left', extent=(-180,180,-180,180),vmin=0.0,vmax=5.0)
            continue
        #Load sequence
        seq=7*nm

        #make fastas for 7mers
        fastafile = run_dir+'poly%s.fasta'%(nm)
        file1=open(fastafile,'w')
        print >>file1,">%s\n%s" % (nm,seq)
        file1.close()

        #run upside_config, nm_hbond1_hbond2_sidechainscale_inversescale_inverseradiusscale
        config = run_dir+'%s_%s.h5'%(nm,prefix)

        GROUPTORUN.append(config)

        if os.path.exists(config): continue   

        upside_config(fastafile, config, 1,
            hbond=False,reference_rama=None,#params_dir+'reference_state_rama.pkl',
            sidechain=None,rotamer=False)

            #edit upside_config
        #z=3
        file2=tb.open_file(config,'r+')
        for z in [1,2,3,4,5]:
            if z==3:
                file2.root.input.potential.rama_map_pot.rama_pot[z] = 0.
                file2.root.input.pivot_moves.proposal_pot[z] = 0.
            else:
                file2.root.input.potential.rama_map_pot.rama_pot[z] = AVERAGERAMA[homorama]
                file2.root.input.pivot_moves.proposal_pot[z] = AVERAGERAMA[homorama]
        file2.close()
        if len(GROUPTORUN)>1:
            #run upside
            while int(os.popen('squeue|grep fillin|wc -l').read().strip())>499:
                time.sleep(30)
                print 'still too many in queue'

            swaps=swap_table2d(1,len(GROUPTORUN))

            print swaps

            nthreads=len(GROUPTORUN)

            jobs[str(GROUPTORUN)] = run_upside('fillin', GROUPTORUN, 120000, 0.2,  n_threads=nthreads,
                              pivot_interval=0.5, replica_interval=1, #time_step=0.007, #max_temp=0.9, 
                              log_level='detailed',estimate_rama_bandwidth=0.5,seed=None,swap_sets=swaps)

            print jobs[str(GROUPTORUN)]
        elif len(GROUPTORUN)==1:
            #run upside
            while int(os.popen('squeue|grep fillin|wc -l').read().strip())>499:
                time.sleep(30)
                print 'still too many in queue'

            jobs[str(GROUPTORUN)] = run_upside('fillin', GROUPTORUN, 120000, 0.2,  n_threads=1,
                              pivot_interval=0.5, #time_step=0.007, #max_temp=0.9, 
                              log_level='detailed',estimate_rama_bandwidth=0.5,seed=None,swap_sets=None)

            print jobs[str(GROUPTORUN)]
