use clap::Parser;
use std::fs::File;
use std::io::Write;
use std::sync::{Arc, Mutex};
use memmap2::Mmap;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std::path::Path;
use std::time::Instant;

#[derive(Parser)]
#[command(name="zip_rfcsv", author="hyf", version="1.0", about="to zip part of RNA framework output file", long_about = None)]
struct Cli {
    #[arg(short,long)]
    input: String,
    #[arg(short,long)]
    output: String,
    #[arg(short,long)]
    thread: usize,
    #[arg(short,long)]
    strand: char,
}

const BAR_LAB: &str = "-\\|/";

fn main() {
    let start = Instant::now();
    let cli = Cli::parse();
    let file = File::open(&cli.input).expect("Cannot open input file please check filename\n");
    let mmap = unsafe { Mmap::map(&file).expect("Cannot mmap file") };
    let mmap = Arc::new(mmap);

    let mut chra:Vec<usize>=Vec::new();
    let mut chrb:Vec<usize>=Vec::new();
    

    println!("Cutting file...");
    let mut regions = Vec::new();
    let mut region_start = 0;
    let mut in_region = false;
    let mut c:usize = 0;
    for (i, &byte) in mmap.iter().enumerate() {
        if !in_region && (byte as char).is_ascii_alphabetic() {
            region_start = i;
            in_region = true;
        } else if in_region && byte == b'\n' {
            if i + 1 < mmap.len() && (mmap[i + 1] as char).is_ascii_alphabetic() {
                regions.push((region_start, i));
                in_region = false;
                if !mmap[region_start..i].contains(&44){
                    chra.push(c+1);
                }
            } else {
                c -= 1;
            }
            c += 1;
        }
    }
    if in_region {
        regions.push((region_start, mmap.len() - 1));
    }
    for i in 0..chra.len()-1 {
        chrb.push(&chra[i+1]-1);
    }
    chrb.push(c+1);
    let mut chunkinfo:Vec<(usize,usize)> = Vec::new();
    for i in 0..chra.len() {
        chunkinfo.push((chra[i],chrb[i]));
    }


    println!("Processing data in parallel...");

    let total_result = Arc::new(Mutex::new(Vec::new()));
    let titleset = Arc::new(Mutex::new(Vec::new()));
    let pool = ThreadPoolBuilder::new().num_threads(cli.thread).build().unwrap();
    pool.install(|| {
        let tchunkinfo = Arc::new(chunkinfo);
        let clen = tchunkinfo.len();
        let tregions = Arc::new(regions);
        let shared_result = Arc::clone(&total_result);
        let shared_titleset = Arc::clone(&titleset);
        let count = Arc::new(Mutex::new(0));
        tchunkinfo.par_iter().for_each(|&(start, end)|{
            let counter = Arc::clone(&count);
            let pside:usize =(start - 1) as usize;
            let qside:usize = end as usize;
            let mut local_result: Vec<String> = Vec::new();
            let mut title = String::new();
            for (index,content) in tregions[pside..qside].iter().enumerate(){
                if index == 0 {
                    title = String::from_utf8_lossy(&mmap[content.0..=content.1]).to_string();
                    let mut titleset = shared_titleset.lock().unwrap();
                    titleset.push(title.clone());
                } else {
                    if String::from_utf8_lossy(&mmap[content.0..=content.1]).to_string().contains(",0,0") {
                        continue;
                    }
                    local_result.push(title.clone().replace('\n', ",")+ &index.clone().to_string() + "," + &String::from_utf8_lossy(&mmap[content.0..=content.1]).to_string().replace('\n', ""));
                }
                
            }
            
            let mut total_result = shared_result.lock().unwrap();
            total_result.push(local_result);
            let mut num = counter.lock().unwrap();
            *num += 1;
 
            print!(
                "\r{}{}{}%",
                "#".repeat(((*num * 100) / clen)/ 2),
                BAR_LAB.chars().nth(((*num * 100) / clen) % 4).unwrap(),
                ((*num * 100) / clen)
            );

        });

    });

    let final_result = total_result.lock().unwrap();
    let final_titleset = titleset.lock().unwrap();

    println!("\nOutput data...");
    //add the file title and output

    let mut file = File::create(cli.output).expect("Unable to create file");
    
    let formatted_titles: String = final_titleset.iter()
        .map(|s| s.replace('\n', ""))
        .collect::<Vec<_>>()
        .join(" ");
    let strand_info: String = final_titleset.iter()
        .map(|_| cli.strand.to_string())
        .collect::<Vec<_>>()
        .join(" ");
    let input_file_name = Path::new(&cli.input).file_name().unwrap().to_str().unwrap();
    let headers = vec![
        "@ColNum 5".to_string(),
        "@ChrID 1".to_string(),
        "@ChrPos 2".to_string(),
        "@Base 3".to_string(),
        format!("@RT_{} 4", input_file_name),
        format!("@BD_{} 5", input_file_name),
        format!("@ChrID_index {}",formatted_titles),
        format!("@ChrID_Strand {}", strand_info),
    ];
    for header in headers {
        writeln!(file, "{}", header).expect("Unable to write header");
    }
    for result in final_result.iter() {
        for line in result {
            writeln!(file, "{}", line).expect("Unable to write data");
        }
    }
    let duration = start.elapsed();
    println!("Total run time: {:?}",duration);
}
