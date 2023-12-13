// Import Modules
include {FAST5_to_POD5; BASECALL; GATHER_BASECALL} from '../modules/basecall'

workflow NANOPORE_STEP_1 {
        
    take:
        basecall_dir
        speed
        modifications

    main:


    BASECALL(basecall_dir, speed, modifications)
    GATHER_BASECALL(id, BASECALL.out.fastq, BASECALL.out.txt)

}
