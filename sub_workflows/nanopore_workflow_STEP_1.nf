// Import Modules

include {FAST5_to_POD5; BASECALL_CPU; BASECALL_CPU_DEMUX} from '../modules/basecall.nf'

workflow NANOPORE_STEP_1 {
        
    take:
        pod5_path
        fast5_path
        speed
        modifications
        config
        trim
        quality_score
        trim_barcode

    main:


    FAST5_to_POD5(fast5_path)
    pod5_path = FAST5_to_POD5.out.mix(pod5_path)


    if (params.basecall_compute?.equalsIgnoreCase("cpu")) {
        
        if (params.basecall_demux == true) {

           println("cpu_demux")
           BASECALL_CPU_DEMUX(pod5_path, speed, modifications, config, trim, quality_score, trim_barcode)
        
        } else {
            
            println("cpu")
            BASECALL_CPU(pod5_path, speed, modifications, config, trim, quality_score)

        }
    
    } else if (params.basecall_compute?.equalsIgnoreCase("gpu")) {
        
        if (params.basecall_demux == true) {

            println("gpu_demux")
            // BASECALL_GPU_DEMUX(pod5_path, speed, modifications, config, trim, quality_score, trim_barcode)

        } else {

            println("gpu")
            BASECALL_GPU(pod5_path, speed, modifications, config, trim, quality_score, trim_barcode)

        }

    }
}
