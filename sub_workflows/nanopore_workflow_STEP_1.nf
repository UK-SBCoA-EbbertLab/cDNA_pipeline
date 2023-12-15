// Import Modules
include {FAST5_to_POD5; BASECALL; GATHER_BASECALL} from '../modules/basecall'

workflow NANOPORE_STEP_1 {
        
    take:
        pod5_path
        fast5_path
        speed
        modifications

    main:


    FAST5_to_POD5(fast5_path)
    pod5_path = FAST5_to_POD5.out.mix(pod5_path)

    BASECALL(pod5_path, speed, modifications)
    // GATHER_BASECALL(BASECALL.out.fastq, BASECALL.out.txt)

}
