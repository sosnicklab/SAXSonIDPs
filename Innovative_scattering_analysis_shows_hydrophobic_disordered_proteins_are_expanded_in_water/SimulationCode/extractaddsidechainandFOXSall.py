from LOADRUNREQs import *


parser = argparse.ArgumentParser(description='Prepare input file',
            usage='use "%(prog)s --help" for more information')
parser.add_argument('--configi', required=True,
            help='config file name without rep.h5')
args = parser.parse_args()
configi=args.configi  
commandargs=[]
for rep in range(0,5):
    print configi
    config='%s_%s_sim.h5' % (configi,rep)
    print config
    if os.path.isdir('%stmpFOXSall/' % config):
        print 'ERROR EXISTS',config
        continue
    commandargs.extend(['cd /fillin/runupside/','&&',upside_dir+'obj/extract_vtf',' --stride 1','--start 50',config,config+'.vtf','&&'
                         ,'vmd',config+'.vtf','<<< \'animate write pdb',config+'.pdb\''
                         ,'&&','mkdir %stmpFOXSall/' % config,'&&','cd %stmpFOXSall/' % config,'&&','csplit %s.pdb /END/ {*}' % config
                         ,'&&',
                         'mkdir fragments','&&','mkdir radii','&&','mkdir rotamers','&&','ln ../fragments/* fragments/','&&','ln ../radii/* radii/','&&','ln ../rotamers/* rotamers/'
,'&&','cp ../SCATD SCATD','&&','cp /fillin/programs/FoXS/foxs foxs'
                         ,'&&',
                        'for f in *; do if [ "${f:0:2}" != \'xx\' ]; then continue; fi; echo "$f"; ./SCATD -i $f -o $f.pdb;'
                        ,'./foxs -s 80 -q 0.4 --excluded_volume 1.0 --water_layer_c2 2.0 $f.pdb;','done','&&','rm -r fragments','&&','rm -r radii','&&','rm -r rotamers','&&'
                        ,'rm SCATD','&&','rm foxs','&&','cd;'])


output_path=config+'_vtfetcFOXSall.output'
args = ['sbatch', '-p', 'amd', '--time=0-%i'%36, '--ntasks=1', 
        '--cpus-per-task=%i'%1, '--export=OMP_NUM_THREADS=%i'%1,
         '--output=%s'%output_path, '--parsable', '--wrap', ' '.join(commandargs)]
job = sp.check_output(args).strip()
