from LOADRUNREQs import *

proteins=[
     ('A','PNt','MP0p7',
      [0., 0.1, 0.399, 0.576, 0.67, 0.751, 0.789, 0.84, 0.872, 0.896, 0.921, 0.943, 0.959, 0.986, 1.003, 1.012, 1.02, 1.033, 1.038, 1.051, 1.062, 1.077, 1.084, 1.088, 1.091, 1.094, 1.096, 1.097, 1.099, 1.1]
      ,[2.38617,0.,0.,0.,0.0080423,-0.0205101,0.0424275,-0.0900227,0.245588,-1.29597,-0.829119,-0.614417,-0.340401,-0.17982,-0.16531,0.0826548,-0.16531]
      ,[334]
     ),
    ('A','polyA','MP0p7',
      [0., 0.1, 0.399, 0.576, 0.67, 0.751, 0.789, 0.84, 0.872, 0.896, 0.921, 0.943, 0.959, 0.986, 1.003, 1.012, 1.02, 1.033, 1.038, 1.051, 1.062, 1.077, 1.084, 1.088, 1.091, 1.094, 1.096, 1.097, 1.099, 1.1]
      ,[2.38617,0.,0.,0.,0.0080423,-0.0205101,0.0424275,-0.0900227,0.245588,-1.29597,-0.829119,-0.614417,-0.340401,-0.17982,-0.16531,0.0826548,-0.16531]
      ,[334]
     ),
    ('A','polyG','MP0p7',
      [0., 0.1, 0.399, 0.576, 0.67, 0.751, 0.789, 0.84, 0.872, 0.896, 0.921, 0.943, 0.959, 0.986, 1.003, 1.012, 1.02, 1.033, 1.038, 1.051, 1.062, 1.077, 1.084, 1.088, 1.091, 1.094, 1.096, 1.097, 1.099, 1.1]
      ,[2.38617,0.,0.,0.,0.0080423,-0.0205101,0.0424275,-0.0900227,0.245588,-1.29597,-0.829119,-0.614417,-0.340401,-0.17982,-0.16531,0.0826548,-0.16531]
      ,[334]
     )]
for nmi,homorama,TYPEi,AMPS,PARAMS,SIZES in proteins:
    DO=True
    for size in [334]:#[50, 73, 107, 156, 228, 334, 488, 714, 1044]:
        for rep in range(0,5):
            config=run_dir+'%s_%s_%s_%s_%s_%s_sim.h5'%(nmi,homorama,size,0,TYPEi,rep)
            try:
                with tb.open_file(config) as t:
                    if len(t.root.output.time[:])!=550:
                        print 'NOT DONE',config
                        DO=False
            except tb.exceptions.NoSuchNodeError:
                print 'NOT STARTED',config
    if os.path.exists(run_dir+'%s_%s_%s_Rgi'%(nmi,homorama,TYPEi)):DO=False
    if DO:
        print 'starting '+run_dir+'%s_%s_%s_Rgi'%(nmi,homorama,TYPEi)
        fileRgi=open(run_dir+'%s_%s_%s_Rgi'%(nmi,homorama,TYPEi),'w')
        fileReei=open(run_dir+'%s_%s_%s_Reei'%(nmi,homorama,TYPEi),'w')
        Rgi=''
        Reei=''
        for numEng in range(0,len(AMPS)):
            Rgi+='\n%s\n' % numEng
            Reei+='\n%s\n' % numEng
            print numEng
            for size in [334]:#[50, 73, 107, 156, 228, 334, 488, 714, 1044]:
                config=run_dir+'%s_%s_%s_%s_%s_%s_sim.h5'%(nmi,homorama,size,numEng,TYPEi,'%s')
                meanRg=[]
                meanRee=[]
                for rep in range(0,5):
                    with tb.open_file(config % rep) as t:
                        try:
                            if ((t.root.output.pos[50:,0]).sum(axis=-1).mean(axis=-1)**2).mean()>0.001:
                                print 'PROBLEM WITH CENTERING!'
                            meanRg.append(np.sqrt((t.root.output.pos[50:,0]**2).sum(axis=-1).mean()))
                            meanRee.append(np.sqrt(((t.root.output.pos[50:,0,-3]-t.root.output.pos[50:,0,0])**2).sum(axis=-1).mean()))
                        except tb.NoSuchNodeError:
                            print 'NOT STARTED YET'
                Rgi+='%s %.5f %.5f\n' % (size,np.array(meanRg).mean(),np.array(meanRg).std()/np.sqrt(len(meanRg)))
                Reei+='%s %.5f %.5f\n' % (size,np.array(meanRee).mean(),np.array(meanRee).std()/np.sqrt(len(meanRee)))
        fileRgi.write(Rgi)
        fileRgi.close()
        fileReei.write(Reei)
        fileReei.close()
    
    for size in SIZES:#[50, 73, 107, 156, 228, 334, 488, 714, 1044]:
        #prefix='%s_%s' % (nmi,rep)
        nm=nmi+'_'+homorama+'_'+str(size)
        if os.path.exists('%s/%s_%s_SCATTERING' % (run_dir,nm,TYPEi)):continue
        commandargs=['python','/fillin/processdataintomathematicaformat.py','--configA=%s/%s' % (run_dir,nm),'--configB=%s' % TYPEi,'--configC=FOXSall','--Params=%s' % len(AMPS)]
        args = ['sbatch', '-p', 'fillin', '--time=0-%i'%36, '--ntasks=1',
        '--cpus-per-task=%i'%1, '--export=OMP_NUM_THREADS=%i'%1,
         '--output=%s'%(run_dir+'%s_%s_%s_LAUNCH'%(nmi,homorama,TYPEi)), '--parsable', '--wrap', ' '.join(commandargs)]
        job = sp.check_output(args).strip()
