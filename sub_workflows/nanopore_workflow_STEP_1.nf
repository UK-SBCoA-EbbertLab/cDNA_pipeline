// Import Modules
include {FAST5_to_POD5; BASECALL; GATHER_BASECALL} from '../modules/basecall'

workflow NANOPORE_STEP_1 {
        
    take:
        basecall_path
        speed
        modifications

    main:


    fast5_dirs = basecall_path.filter{tuple -> tuple[1].toString().contains('/fast5_pass/')}
    pod5_dirs = basecall_path.filter{tuple -> tuple[1].toString().contains('/pod5_pass/')}
    

    FAST5_to_POD5(fast5_dirs)
    basecall_path = FAST5_to_POD5.out.mix(pod5_dirs)

    BASECALL(basecall_path, speed, modifications)
    GATHER_BASECALL(BASECALL.out.fastq, BASECALL.out.txt)

}
