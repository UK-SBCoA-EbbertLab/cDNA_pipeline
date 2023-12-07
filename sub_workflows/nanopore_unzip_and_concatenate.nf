// Import Modules
include {UNZIP_AND_CONCATENATE} from '../modules/unzip_and_concatenate.sh'

workflow NANOPORE_UNZIP_AND_CONCATENATE {
        
    take:
        path

    main:


    UNZIP_AND_CONCATENATE(path)

}
