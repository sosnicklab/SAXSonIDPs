from LOADRUNREQs import *


proti=1
#MAke Average PNt RAMAs
#jobs = dict()
AVERAGERAMA={}
RAMATYPE={}
proteins=[('PNt','AverageRama'),('polyG','AverageRama'),('polyA','AverageRama')]

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
                    tmpRAMA=AVERAGERAMA[nm]
                    
                    continue


                print 'ERROR DOES NOT EXIST RUN fileA',nm,TYPE,PPII,LJDFIRE


#GOTOTHISPART

proteins=[
('A','PNt','MP0p7',
      [0., 0.1, 0.399, 0.576, 0.67, 0.751, 0.789, 0.84, 0.872, 0.896, 0.921, 0.943, 0.959, 0.986, 1.003, 1.012, 1.02, 1.033, 1.038, 1.051, 1.062, 1.077, 1.084, 1.088, 1.091, 1.094, 
1.096, 1.097, 1.099, 1.1]
      ,[2.38617,0.,0.,0.,0.0080423,-0.0205101,0.0424275,-0.0900227,0.245588,-1.29597,-0.829119,-0.614417,-0.340401,-0.17982,-0.16531,0.0826548,-0.16531]
),
('A','polyA','MP0p7',
      [0., 0.1, 0.399, 0.576, 0.67, 0.751, 0.789, 0.84, 0.872, 0.896, 0.921, 0.943, 0.959, 0.986, 1.003, 1.012, 1.02, 1.033, 1.038, 1.051, 1.062, 1.077, 1.084, 1.088, 1.091, 1.094, 
1.096, 1.097, 1.099, 1.1]
      ,[2.38617,0.,0.,0.,0.0080423,-0.0205101,0.0424275,-0.0900227,0.245588,-1.29597,-0.829119,-0.614417,-0.340401,-0.17982,-0.16531,0.0826548,-0.16531]
),
('A','polyG','MP0p7',
      [0., 0.1, 0.399, 0.576, 0.67, 0.751, 0.789, 0.84, 0.872, 0.896, 0.921, 0.943, 0.959, 0.986, 1.003, 1.012, 1.02, 1.033, 1.038, 1.051, 1.062, 1.077, 1.084, 1.088, 1.091, 1.094, 
1.096, 1.097, 1.099, 1.1]
      ,[2.38617,0.,0.,0.,0.0080423,-0.0205101,0.0424275,-0.0900227,0.245588,-1.29597,-0.829119,-0.614417,-0.340401,-0.17982,-0.16531,0.0826548,-0.16531]
)]
cnt=0
for nmi,homorama,TYPEi,AMPS,PARAMS in proteins:
    #clear_output()
    ALLreps={}
    for rep in range(0,5):
        TYPE='%s_%s' % (TYPEi,rep)
        if RAMATYPE[homorama]=='AverageRamaB':TYPE='%sb_%s' % (TYPEi,rep)
        for size in [334]:
            clear_output()
            ALLreps.setdefault(size,{})
            GROUPTORUN=[]
            run=True
            cnt+=1
            if cnt>2:
                clear_output()
                cnt=0
                exit()
            for MPDFIRE in range(0,len(AMPS)):
                #prefix='%s_%s' % (nmi,rep)
                ALLreps[size].setdefault(MPDFIRE,[])
                nm=nmi+'_'+homorama+'_'+str(size)

                config=run_dir+'%s_%s_%s_sim.h5'%(nm,MPDFIRE,TYPE)

                #ramafile=run_dir+'%s_%s.rama.pkl'%(nmi,TYPE)
                if os.path.exists(config):
                    #SizeandRg(config,size)
                    if size==49:
                        if len(ALLreps[size][0])==0:
                            try :
                                print config
                                #RgandDist(config,50)
                                #Temp(config,50)
                            except tb.HDF5ExtError:
                                print 'OpenIssue'
                            print
                    with tb.open_file(config) as t:
                        try:
                            ALLreps[size][MPDFIRE].append(np.sqrt((t.root.output.pos[50::,0]**2).sum(axis=-1).mean(axis=-2).mean()))
                        except tb.NoSuchNodeError:
                            print 'NOT STARTED YET'
                    continue
                    #run=False
                    continue

                #Load sequence
                seq=size*nmi

                #make fastas for n-mers
                fastafile = run_dir+'%s.fasta'%(nm)
                file1=open(fastafile,'w')
                print >>file1,">%s\n%s" % (nm,seq)
                file1.close()

                print 'STARTING FINAL',nm,TYPE,MPDFIRE

                sidechaini='Cbhomonu_0.h5'
                if MPDFIRE==None:
                    sidechaini=None

                upside_config(run_dir+'%s.fasta' % nm, config, 1,
                      hbond=False, reference_rama=None, 
                      sidechain=sidechaini,rama_pot='rama_libraries.h5',rama_pot_product=True,#sidechain_scale=MJi, inverse_scale=opts[1], inverse_radius_scale=opts[2],
                      rotamer=False)

                #make rama

                print 'making rama',nm,TYPE,MPDFIRE

                #Load TrueRama

                file2=tb.open_file(config,'r+')

                #Load sequence

                ALLRAMA=np.ndarray(shape=(size,72,72))

                ramax = run_dir+'%s_%s_%s.h5.rama.pkl'%(nmi,homorama,rep)
                if RAMATYPE[homorama]=='AverageRamaB':
                    ramax = run_dir+'%s_%s_%sb.h5.rama.pkl'%(nmi,homorama,rep)
                
                ALLRAMAi=np.swapaxes(cp.load(open(ramax))[3],0,1)
                for x in range(1,size-1):
                    #ALLRAMA[x]=ALLRAMAi
                    file2.root.input.potential.rama_map_pot.rama_pot[x] =AVERAGERAMA[homorama]+np.log(ALLRAMAi/sum(sum(ALLRAMAi))+0.000000001)
                    file2.root.input.pivot_moves.proposal_pot[x] = AVERAGERAMA[homorama]+np.log(ALLRAMAi/sum(sum(ALLRAMAi))+0.000000001)
                
                file2.root.input.potential.radial.interaction_param[:,:]=[PARAMS[0]]+list(AMPS[MPDFIRE]*np.array(PARAMS[1:]))

                
                file2.close()

                GROUPTORUN.append(config)

            if run==False:continue
            if len(GROUPTORUN)<1:continue
            print nm,len(GROUPTORUN)
            swaps=swap_table2d(1,len(GROUPTORUN))

            print swaps

            nthreads=(len(GROUPTORUN)+2)/3

            print nthreads
            if len(GROUPTORUN)>1:
                jobs[str(GROUPTORUN)] = run_upside('fillin', GROUPTORUN, 11000, 20,  n_threads=nthreads,
                                  pivot_interval=0.5, replica_interval=1,# #time_step=0.007, #max_temp=0.9, 
                                  log_level='detailed',seed=None ,estimate_rama_bandwidth=None ,swap_sets=swaps)
            else :
                jobs[str(GROUPTORUN)] = run_upside(
                                    'fillin'
                                    , GROUPTORUN, 11000, 20,  n_threads=nthreads,
                                  pivot_interval=0.5, #replica_interval=1,# #time_step=0.007, #max_temp=0.9, 
                                  log_level='detailed',seed=None ,estimate_rama_bandwidth=None) #,swap_sets=swaps)

            print jobs[str(GROUPTORUN)]
    #print 'index','amp','Rg156','Rg107','nu_aprox'
    KEYS=ALLreps.keys()
    MDFIRES=ALLreps[KEYS[0]].keys()
    for x in MDFIRES:
        TMPSTRING='%s  ' % AMPS[x]
        for y in [334]:#[50, 73, 107, 156, 228, 334, 488, 714, 1044]:
            if not y in KEYS:continue
            TMPSTRING+='%.2f  ' % np.array(ALLreps[y][x]).mean()
        print TMPSTRING#AMPS[z],round(FIRST[z],1)
        #print z,AMPS[z],round(SECOND[z],1),round(FIRST[z],1),round(np.log(SECOND[z]/FIRST[z])/np.log(156./107.),2)
            #os.remove(ramafile)
