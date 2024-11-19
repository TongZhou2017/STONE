use std::collections::HashMap;
use std::fs::File;
use clap::Parser;
use std::{char, str};
use std::sync::{Arc, Mutex};
use std::io::{BufReader, BufRead, Write};
use std::time::Instant;
use std::collections::HashSet;
use rand::Rng;
use rayon::ThreadPoolBuilder;
use rayon::iter::ParallelIterator;
use rayon::iter::IntoParallelRefIterator;

fn chr_to_nc(chr: &str) -> String {
    let mapping = HashMap::from([
        ("chr1", "NC_000001"),
        ("chr2", "NC_000002"),
        ("chr3", "NC_000003"),
        ("chr4", "NC_000004"),
        ("chr5", "NC_000005"),
        ("chr6", "NC_000006"),
        ("chr7", "NC_000007"),
        ("chr8", "NC_000008"),
        ("chr9", "NC_000009"),
        ("chr10", "NC_000010"),
        ("chr11", "NC_000011"),
        ("chr12", "NC_000012"),
        ("chr13", "NC_000013"),
        ("chr14", "NC_000014"),
        ("chr15", "NC_000015"),
        ("chr16", "NC_000016"),
        ("chr17", "NC_000017"),
        ("chr18", "NC_000018"),
        ("chr19", "NC_000019"),
        ("chr20", "NC_000020"),
        ("chr21", "NC_000021"),
        ("chr22", "NC_000022"),
        ("chrX", "NC_000023"),
        ("chrY", "NC_000024"),
        ("chrM", "NC_012920"),  
    ]);

    mapping.get(chr).map(|s| s.to_string()).unwrap_or_else(|| chr.to_string())
}

fn chr_to_nc_mouse(chr: &str) -> String {
    let mapping = HashMap::from([
        ("chr1", "NC_000067"),
        ("chr2", "NC_000068"),
        ("chr3", "NC_000069"),
        ("chr4", "NC_000070"),
        ("chr5", "NC_000071"),
        ("chr6", "NC_000072"),
        ("chr7", "NC_000073"),
        ("chr8", "NC_000074"),
        ("chr9", "NC_000075"),
        ("chr10", "NC_000076"),
        ("chr11", "NC_000077"),
        ("chr12", "NC_000078"),
        ("chr13", "NC_000079"),
        ("chr14", "NC_000080"),
        ("chr15", "NC_000081"),
        ("chr16", "NC_000082"),
        ("chr17", "NC_000083"),
        ("chr18", "NC_000084"),
        ("chr19", "NC_000085"),
        ("chrX", "NC_000086"),
        ("chrY", "NC_000087"),
        ("chrM", "NC_005089"),  
    ]);

    mapping.get(chr).map(|s| s.to_string()).unwrap_or_else(|| chr.to_string())
}

fn chr_to_nc_yeast(chr: &str) -> String {
    let mapping = HashMap::from([
        ("chrI", "NC_001133"),
        ("chrII", "NC_001134"),
        ("chrIII", "NC_001135"),
        ("chrIV", "NC_001136"),
        ("chrV", "NC_001137"),
        ("chrVI", "NC_001138"),
        ("chrVII", "NC_001139"),
        ("chrVIII", "NC_001140"),
        ("chrIX", "NC_001141"),
        ("chrX", "NC_001142"),
        ("chrXI", "NC_001143"),
        ("chrXII", "NC_001144"),
        ("chrXIII", "NC_001145"),
        ("chrXIV", "NC_001146"),
        ("chrXV", "NC_001147"),
        ("chrXVI", "NC_001148"),
        ("chrM", "NC_001224"), 
    ]);

    mapping.get(chr).map(|s| s.to_string()).unwrap_or_else(|| chr.to_string())
}

fn chr_to_nc_at(chr: &str) -> String {
    let mapping = HashMap::from([
        ("chr1", "NC_003070"),
        ("chr2", "NC_003071"),
        ("chr3", "NC_003074"),
        ("chr4", "NC_003075"),
        ("chr5", "NC_003076"),
        ("chrM", "NC_001284"),  
        ("chrC", "NC_000932"),  
    ]);

    mapping.get(chr).map(|s| s.to_string()).unwrap_or_else(|| chr.to_string())
}

fn get_key_prefix(key: &str) -> &str {
    key.split('.').next().unwrap_or(key)
}

#[derive(Debug)]
struct Record {
    strand: String,
    position: usize,
    info: String,  
}

#[derive(Debug)]
struct BedRecord {
    start: usize,
    end: usize,
    strand: String,
    name: String,
    tname: String,
}

fn get_chr_nc_for_species(species: &str, chr: &str) -> String {
    match species {
        "hu" => chr_to_nc(chr),
        "mo" => chr_to_nc_mouse(chr),
        "ye" => chr_to_nc_yeast(chr),
        "at" => chr_to_nc_at(chr),
        _ => chr.to_string(),  
    }
}

fn read_bed_file(file_path: &str,species: &str) -> HashMap<String, Vec<BedRecord>> {
    let file = File::open(file_path).unwrap();
    let reader = BufReader::new(file);

    let mut records_map: HashMap<String, Vec<BedRecord>> = HashMap::new();

    for line in reader.lines() {
        let line = line.unwrap();
        let fields: Vec<&str> = line.split_whitespace().collect();

        if fields.len() >= 1 {
            let chr = fields[0].to_string();                // 
            let start: usize = fields[1].parse().unwrap();  
            let end: usize = fields[2].parse().unwrap();    // 
            let strand = fields[3].to_string();             //             
            let name = fields[4].split('=').next().unwrap_or(fields[4]).to_string(); 
            let tname = fields[5].to_string();


            let nc_chr = get_chr_nc_for_species(&species, &chr);

            let record = BedRecord {
                start,
                end,
                strand,
                name,
                tname,
            };

            records_map.entry(nc_chr).or_insert_with(Vec::new).push(record);
        }
    }

    records_map
}


fn read_txt(file_path: &str) -> HashMap<String, Vec<Record>> {
    let file = File::open(file_path).unwrap();
    let reader = BufReader::new(file);

    let mut records_map: HashMap<String, Vec<Record>> = HashMap::new();

    for (i, line) in reader.lines().enumerate() {
        let line = line.unwrap();
        if i == 0 { continue; } 

        let fields: Vec<&str> = line.split(',').collect();

        let chrid = fields[0].to_string();
        let strand = fields[1].to_string();
        let position = fields[2].parse::<usize>().unwrap();
        let info = fields[3..].join(","); 

        let rec = Record {
            strand,
            position,
            info,
        };

        records_map.entry(chrid).or_insert_with(Vec::new).push(rec);
    }

    records_map
}

#[derive(Parser)]
#[command(name="bgsg", author="hyf", version="1.0", about="bgsg species:\nhuman-hu\nmouse-mo\nyeast-ye\nArabidopsis thaliana-at", long_about = None)]
struct Cli {
    #[arg(short,long)]
    mergepath: String,
    #[arg(short,long)]
    bedpath: String,
    #[arg(short,long)]
    output: String,
    #[arg(short,long)]
    species: String,
}


fn main() {
    let start = Instant::now();
    let cli = Cli::parse();
    let species = &cli.species;
    println!("load file");
    let mergeset = Arc::new(read_txt(&cli.mergepath));
    let bedset = Arc::new(read_bed_file(&cli.bedpath,species));
    let file = Arc::new(Mutex::new(File::create(cli.output).expect("Unable to create file")));
    
    {
        let mut file_guard = file.lock().unwrap();
        writeln!(file_guard, "ChrID,geneid,transcriptid,position,pipe_truncation_Strand,rf_mutation_Base,rf_mutation_Count,rf_mutation_Depth,rf_mutation_AC,rf_mutation_AG,rf_mutation_AT,rf_mutation_CA,rf_mutation_CG,rf_mutation_CT,rf_mutation_GA,rf_mutation_GC,rf_mutation_GT,rf_mutation_TA,rf_mutation_TC,rf_mutation_TG,rf_mutation_ins,rf_mutation_del,pipe_truncation_Base,pipe_truncation_count,pipe_truncation_BD,base_A,base_C,base_G,base_T,modified_string").expect("Unable to write header");
    }

    let thread_count = 32;

    println!("parallel");
    let pool = ThreadPoolBuilder::new().num_threads(thread_count).build().unwrap();
    pool.install(|| {
        bedset.par_iter().for_each(|(bed_key, bed_records)| {
            let bed_key_prefix = get_key_prefix(bed_key);
            mergeset.iter().for_each(|(merge_key, records)| {
                let merge_key_prefix = get_key_prefix(merge_key);

                if bed_key_prefix == merge_key_prefix {
                    let mut local_result = String::new();
                    process_records(bed_key_prefix, bed_records, records, &mut local_result);


                    let mut file_guard = file.lock().unwrap();
                    writeln!(file_guard, "{}", local_result.trim()).expect("Unable to write data");

                    local_result.clear(); 
                }
            });
        });
    });

    let duration = start.elapsed();
    println!("gene time elapsed: {:?}", duration);
}


fn process_records(bedkeyprefix: &str, bed_records: &Vec<BedRecord>, records: &Vec<Record>, result: &mut String) {
    let mut rng = rand::thread_rng();
    for bed_record in bed_records {
        let start = bed_record.start;
        let end = bed_record.end;

        let mut filtered_records: Vec<&Record> = records.iter()
            .filter(|r| r.strand == bed_record.strand && r.position >= start && r.position <= end)
            .collect();
        let mut seen = HashSet::new();
        filtered_records.retain(|r| seen.insert((&r.strand, r.position)));
        for i in filtered_records {
            let rnum:u8= rng.gen_range(0..=1);
            let chars_as_string: String = parse_and_extend_info(&i.info);
            result.push_str(&format!("{},{},{},{},{},{},{}\n", bedkeyprefix, bed_record.name, bed_record.tname, i.position, bed_record.strand, chars_as_string, rnum));
        }
        

    }
}

fn parse_and_extend_info(info: &str) -> String {
    let fields: Vec<&str> = info.split(',').collect();

    // 提取ref碱基和测序深度及突变数
    let ref_base:char = fields[0].chars().next().unwrap();
    let mutation_count: usize = fields[1].parse().unwrap_or(0);
    let total_depth: usize = fields[2].parse().unwrap_or(0);

    let ref_count = total_depth.saturating_sub(mutation_count);

    let mut a_count = 0;
    let mut c_count = 0;
    let mut g_count = 0;
    let mut t_count = 0;


    let a_indices = vec![3, 4, 5]; 
    let c_indices = vec![6, 7, 8]; 
    let g_indices = vec![9, 10, 11]; 
    let t_indices = vec![12, 13, 14]; 

    match ref_base {
        'A' => {
            a_count = ref_count;
            c_count = fields[a_indices[0]].parse().unwrap();
            g_count = fields[a_indices[1]].parse().unwrap();
            t_count = fields[a_indices[2]].parse().unwrap();
        },
        'C' => {
            c_count = ref_count;
            a_count = fields[c_indices[0]].parse().unwrap();
            g_count = fields[c_indices[1]].parse().unwrap();
            t_count = fields[c_indices[2]].parse().unwrap();
        },
        'G' => {
            g_count = ref_count;
            a_count = fields[g_indices[0]].parse().unwrap();
            c_count = fields[g_indices[1]].parse().unwrap();
            t_count = fields[g_indices[2]].parse().unwrap();
        },
        'T' => {
            t_count = ref_count;
            a_count = fields[t_indices[0]].parse().unwrap();
            c_count = fields[t_indices[1]].parse().unwrap();
            g_count = fields[t_indices[2]].parse().unwrap();
        },
        _ => {},
    }


    let extended_info = format!(
        "{},{},{},{},{}",
        info,
        a_count,
        c_count,
        g_count,
        t_count
    );

    extended_info
}

