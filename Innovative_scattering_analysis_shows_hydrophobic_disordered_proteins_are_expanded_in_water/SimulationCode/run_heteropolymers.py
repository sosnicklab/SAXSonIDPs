from LOADRUNREQs import *

proteins=[
    ('PNt','ALLH','MP0p7',
      [0., 0.423, 0.556, 0.662, 0.731, 0.784, 0.829, 0.857, 0.879, 0.894, 
       0.908, 0.92, 0.941, 0.957, 0.97, 0.979, 0.985, 0.992, 1.004, 1.019, 
       1.028, 1.032, 1.038, 1.046, 1.052, 1.058, 1.065, 1.071, 1.077, 1.082]
      ,[2.38617,0.,0.,0.,0.0080423,-0.0205101,0.0424275,-0.0900227,0.245588,-1.29597,-0.829119,-0.614417,-0.340401,-0.17982,-0.16531,0.0826548,-0.16531]
     ),
    ('FhuA','ALLH','MP0p7',
      [0., 0.423, 0.556, 0.662, 0.731, 0.784, 0.829, 0.857, 0.879, 0.894, 
       0.908, 0.92, 0.941, 0.957, 0.97, 0.979, 0.985, 0.992, 1.004, 1.019, 
       1.028, 1.032, 1.038, 1.046, 1.052, 1.058, 1.065, 1.071, 1.077, 1.082]
      ,[2.38617,0.,0.,0.,0.0080423,-0.0205101,0.0424275,-0.0900227,0.245588,-1.29597,-0.829119,-0.614417,-0.340401,-0.17982,-0.16531,0.0826548,-0.16531]
     ),
     ('PNt','HP','MP0p7',
      list((3.0*i for i in [0., 0.423, 0.556, 0.662, 0.731, 0.784, 0.829, 0.857, 0.879, 0.894, 
       0.908, 0.92, 0.941, 0.957, 0.97, 0.979, 0.985, 0.992, 1.004, 1.019, 
       1.028, 1.032, 1.038, 1.046, 1.052, 1.058, 1.065, 1.071, 1.077, 1.082]))
      ,[2.38617,0.,0.,0.,0.0080423,-0.0205101,0.0424275,-0.0900227,0.245588,-1.29597,-0.829119,-0.614417,-0.340401,-0.17982,-0.16531,0.0826548,-0.16531]
     )]
cnt=0
for nmi,hydroptype,TYPEi,AMPS,PARAMS in proteins:
    #clear_output()
    ALLreps={}
    for rep in range(0,5):
        TYPE='%s_%s' % (TYPEi,rep)
        GROUPTORUN=[]
        run=True
        cnt+=1
        if cnt>2:
            cnt=0
        for MPDFIRE in range(0,len(AMPS)):
            prefix='%s_%s' % (nmi,rep)
            ALLreps.setdefault(MPDFIRE,[])
            nm=nmi+'_'+hydroptype

            config=run_dir+'%s_%s_%s_sim.h5'%(nm,MPDFIRE,TYPE)

            #ramafile=run_dir+'%s_%s.rama.pkl'%(nmi,TYPE)
            if os.path.exists(config):
                #SizeandRg(config,size)
                if len(ALLreps[0])==0:
                    try :
                        print config
                        #RgandDist(config,50)
                        #Temp(config,50)
                    except tb.HDF5ExtError:
                        print 'OpenIssue'
                    print
                with tb.open_file(config) as t:
                    try:
                        ALLreps[MPDFIRE].append(np.sqrt((t.root.output.pos[50::,0]**2).sum(axis=-1).mean(axis=-2).mean()))
                    except tb.NoSuchNodeError:
                        #print t.root.output
                        print 'NOT STARTED YET'
                continue
                #run=False
                continue


            #define fasta file
            fastafile = '%s.fasta'%(nmi)

            print 'STARTING FINAL',nm,TYPE,MPDFIRE

            sidechaini='Cbhomonu_0.h5'
            if MPDFIRE==None:
                sidechaini=None
                
            coillib='rama_libraries.h5'
            if TYPEi=='MP0p7Coil':
                coillib='rama_librariesConly.h5'

            upside_config('%s.fasta' % nmi, config, 1,
                  hbond=False, reference_rama=params_dir+'reference_state_rama.pkl', 
                  sidechain=sidechaini,rama_pot=coillib,rama_pot_product=True,#sidechain_scale=MJi, inverse_scale=opts[1], inverse_radius_scale=opts[2],
                  rotamer=False)

            #make rama

            print 'making rama',nm,TYPE,MPDFIRE

            #Load TrueRama

            file2=tb.open_file(config,'r+')

            #Load sequence

            file2.root.input.potential.radial.interaction_param[:,:]=[PARAMS[0]]+list(AMPS[MPDFIRE]*np.array(PARAMS[1:]))
            
            #Change radial type
            if hydroptype=='HP':
                file2.root.input.potential.radial.interaction_param[0,:]=[PARAMS[0]]+list(0.0*np.array(PARAMS[1:]))
                file2.root.input.potential.radial.interaction_param[:,0]=[PARAMS[0]]+list(0.0*np.array(PARAMS[1:]))
                for x in range(0,len(file2.root.input.potential.rama_map_pot.rama_pot[:])):
                    if file2.root.input.sequence[x] in ['ALA','VAL','LEU','ILE','MET','PHE','TYR','TRP']:
                        file2.root.input.potential.radial.type[x]=1
                    else :
                        file2.root.input.potential.radial.type[x]=0
            elif hydroptype=='4CLUST21':
                file2.root.input.potential.radial.interaction_param[:,:]=[PARAMS[0]]+list(0.0*np.array(PARAMS[1:]))
                file2.root.input.potential.radial.interaction_param[1,1]=[PARAMS[0]]+list(2.*np.array(PARAMS[1:]))
                file2.root.input.potential.radial.interaction_param[2,2]=[PARAMS[0]]+list(2.*np.array(PARAMS[1:]))
                file2.root.input.potential.radial.interaction_param[3,3]=[PARAMS[0]]+list(2.*np.array(PARAMS[1:]))
                file2.root.input.potential.radial.interaction_param[4,4]=[PARAMS[0]]+list(2.*np.array(PARAMS[1:]))
                LENGTHS=[21,62,21,63,21,62,21,63]
                TYPE=[1,0,2,0,3,0,4,0]
                countlength=0
                counttype=0
                for x in range(0,len(file2.root.input.potential.rama_map_pot.rama_pot[:])):
                    countlength+=1
                    if countlength>LENGTHS[counttype]:
                        countlength=1
                        counttype+=1
                    file2.root.input.potential.radial.type[x]=TYPE[counttype]
            elif hydroptype=='4CLUST42':
                file2.root.input.potential.radial.interaction_param[:,:]=[PARAMS[0]]+list(0.0*np.array(PARAMS[1:]))
                file2.root.input.potential.radial.interaction_param[1,1]=[PARAMS[0]]+list(2.*np.array(PARAMS[1:]))
                file2.root.input.potential.radial.interaction_param[2,2]=[PARAMS[0]]+list(2.*np.array(PARAMS[1:]))
                file2.root.input.potential.radial.interaction_param[3,3]=[PARAMS[0]]+list(2.*np.array(PARAMS[1:]))
                file2.root.input.potential.radial.interaction_param[4,4]=[PARAMS[0]]+list(2.*np.array(PARAMS[1:]))
                LENGTHS=[42 , 41 , 42 , 42 , 42 , 41 , 42 , 42]
                TYPE=[1,0,2,0,3,0,4,0]
                countlength=0
                counttype=0
                for x in range(0,len(file2.root.input.potential.rama_map_pot.rama_pot[:])):
                    countlength+=1
                    if countlength>LENGTHS[counttype]:
                        countlength=1
                        counttype+=1
                    file2.root.input.potential.radial.type[x]=TYPE[counttype]
            elif hydroptype=='2CLUST42':
                file2.root.input.potential.radial.interaction_param[:,:]=[PARAMS[0]]+list(0.0*np.array(PARAMS[1:]))
                file2.root.input.potential.radial.interaction_param[1,1]=[PARAMS[0]]+list(2.*np.array(PARAMS[1:]))
                file2.root.input.potential.radial.interaction_param[2,2]=[PARAMS[0]]+list(2.*np.array(PARAMS[1:]))
                LENGTHS=[42 , 125 , 42 , 125]
                TYPE=[1,0,2,0]
                countlength=0
                counttype=0
                for x in range(0,len(file2.root.input.potential.rama_map_pot.rama_pot[:])):
                    countlength+=1
                    if countlength>LENGTHS[counttype]:
                        countlength=1
                        counttype+=1
                    file2.root.input.potential.radial.type[x]=TYPE[counttype]
            elif hydroptype=='2CLUST84':
                file2.root.input.potential.radial.interaction_param[:,:]=[PARAMS[0]]+list(0.0*np.array(PARAMS[1:]))
                file2.root.input.potential.radial.interaction_param[1,1]=[PARAMS[0]]+list(2.*np.array(PARAMS[1:]))
                file2.root.input.potential.radial.interaction_param[2,2]=[PARAMS[0]]+list(2.*np.array(PARAMS[1:]))
                LENGTHS=[84 , 83 , 84 , 83]
                TYPE=[1,0,2,0]
                countlength=0
                counttype=0
                for x in range(0,len(file2.root.input.potential.rama_map_pot.rama_pot[:])):
                    countlength+=1
                    if countlength>LENGTHS[counttype]:
                        countlength=1
                        counttype+=1
                    file2.root.input.potential.radial.type[x]=TYPE[counttype]
            elif hydroptype=='1CLUST84':
                file2.root.input.potential.radial.interaction_param[:,:]=[PARAMS[0]]+list(0.0*np.array(PARAMS[1:]))
                file2.root.input.potential.radial.interaction_param[1,1]=[PARAMS[0]]+list(2.*np.array(PARAMS[1:]))
                LENGTHS=[84 , 250]
                TYPE=[1,0]
                countlength=0
                counttype=0
                for x in range(0,len(file2.root.input.potential.rama_map_pot.rama_pot[:])):
                    countlength+=1
                    if countlength>LENGTHS[counttype]:
                        countlength=1
                        counttype+=1
                    file2.root.input.potential.radial.type[x]=TYPE[counttype]
            elif hydroptype=='1CLUST167':
                file2.root.input.potential.radial.interaction_param[:,:]=[PARAMS[0]]+list(0.0*np.array(PARAMS[1:]))
                file2.root.input.potential.radial.interaction_param[1,1]=[PARAMS[0]]+list(2.*np.array(PARAMS[1:]))
                LENGTHS=[167 , 167]
                TYPE=[1,0]
                countlength=0
                counttype=0
                for x in range(0,len(file2.root.input.potential.rama_map_pot.rama_pot[:])):
                    countlength+=1
                    if countlength>LENGTHS[counttype]:
                        countlength=1
                        counttype+=1
                    file2.root.input.potential.radial.type[x]=TYPE[counttype]
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
    MDFIRES=ALLreps.keys()
    for x in MDFIRES:
        TMPSTRING='%s  ' % AMPS[x]
        TMPSTRING+='%.2f  ' % np.array(ALLreps[x]).mean()
        print TMPSTRING#AMPS[z],round(FIRST[z],1)
        #print z,AMPS[z],round(SECOND[z],1),round(FIRST[z],1),round(np.log(SECOND[z]/FIRST[z])/np.log(156./107.),2)
            #os.remove(ramafile)
