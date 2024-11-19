The work is done in rust, only the source code is provided here, 
to get the binary please install cargo locally and use it

Note that you add the relevant crate to the cargo.toml file!

1.source workflow
                
RNAframework
    (1)csv->[zip_rfcsv]->zippedcsv
    (2)txt->[zip_rftxt]->[zip_rftxt]->zippedtxt

icSHAPE-pipe

    (1)output->[zip_pipe]->zippedpipe

2.Merge
    use the files in source above to generate mergefile

    zippedcsv----|
    zippedpipe-|-> [merge] -> mergedfile ->[mbreport]->statistic results
    zippedtxt----|

3.extract or search
    exact all the genes in bed.file or search single gene through bed.file

               
    mergedfile---> [bgsg] ----> all gene's result


3.Usage

 (1)zip_pipe

    Usage: {} <input_file> <output_file>

 (2)zip_rfcsv:

    to zip part of RNA framework output csv
    Usage: zip_rfcsv --input <INPUT> --output <OUTPUT> --thread <THREAD> --strand <STRAND>
    Options:
    -i, --input <INPUT>    
    -o, --output <OUTPUT>  
    -t, --thread <THREAD>  
    -s, --strand <STRAND>  
    -h, --help             Print help
    -V, --version          Print version

 (3)zip_rftxt

    Usage: {} <input_file> <output_file> <num_threads>

 (4)zip_rftxt2
    
    to reduce previous step of RNA framework output file
    Usage: zip_rftxt --input <INPUT> --output <OUTPUT> --thread <THREAD> --strand <STRAND>
    Options:
    -i, --input <INPUT>    
    -o, --output <OUTPUT>  
    -t, --thread <THREAD>  
    -s, --strand <STRAND>  
    -h, --help             Print help
    -V, --version          Print version

(5)Merge three files

    merge three files
    Usage: merge --csv <CSV> --pipe <PIPE> --txt <TXT> --output <OUTPUT>
    Options:
    -c, --csv <CSV>        
    -p, --pipe <PIPE>      
    -t, --txt <TXT>        
    -o, --output <OUTPUT>  
    -h, --help             Print help
    -V, --version          Print version


(6)bgsg

    exact all genes from bedfile
    Usage: bgsg --mergepath <MERGEPATH> --bedpath <BEDPATH> --output <OUTPUT>
    Options:
    -m, --mergepath <MERGEPATH>  
    -b, --bedpath <BEDPATH>      
    -o, --output <OUTPUT>        
    -h, --help                   Print help
    -V, --version                Print version

(5)mbreport

    statistic
    Usage: mbreport --input <INPUT> --depth <DEPTH> --output <OUTPUT>
    Options:
    -i, --input <INPUT>    
    -d, --depth <DEPTH>    #Use this to determine the minimum depth
    -o, --output <OUTPUT>  
    -h, --help             Print help
    -V, --version          Print version


