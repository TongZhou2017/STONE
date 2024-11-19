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
#[command(name="zip_rftxt2", author="hyf", version="1.0", about="to reduction previous step of RNA framework output file", long_about = None)]
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

#[allow(unused_assignments)]
fn main() {
    let start = Instant::now();
    let cli = Cli::parse();
    println!("loading file...");
    let file = File::open(&cli.input).expect("Cannot open input file");
    let mmap = unsafe { Mmap::map(&file).expect("Cannot mmap file") };
    let mmap = Arc::new(mmap);

    let mut chra:Vec<usize>=Vec::new();
    let mut chrb:Vec<usize>=Vec::new();

    println!("Cutting file...");
    let mut regions = Vec::new();
    let mut region_start = 0;
    let mut in_region = false;
    let mut c = 0;
    for (i, &byte) in mmap.iter().enumerate() {
        if !in_region && (byte as char).is_ascii_alphabetic() {
            region_start = i;
            in_region = true;
        } else if in_region && byte == b'\n' {
            if i + 1 < mmap.len() && (mmap[i + 1] as char).is_ascii_alphabetic() {
                regions.push((region_start, i));
                in_region = false;
                if !mmap[region_start..i].contains(&9){
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
        let count = Arc::new(Mutex::new(0));
        let tchunkinfo = Arc::new(chunkinfo);
        let clen = tchunkinfo.len();
        let tregions = Arc::new(regions);
        let shared_result = Arc::clone(&total_result);
        let shared_titleset = Arc::clone(&titleset);
        tchunkinfo.par_iter().for_each(|&(start, end)|{
            let counter = Arc::clone(&count);
            let pside:usize =(start - 1) as usize;
            let qside:usize = end as usize;
            let mut outkv:Vec<Vec<(i32,(usize,usize))>> = Vec::new();
            let mut kv :Vec<(i32,(usize,usize))> = Vec::new();
            let mut title = String::new();
            for (index,content) in tregions[pside..qside].iter().enumerate(){
                if index == 0 {
                    title = String::from_utf8_lossy(&mmap[content.0..=content.1]).to_string();
                    let mut titleset = shared_titleset.lock().unwrap();
                    titleset.push(title.clone());
                } else {
                    println!("{:?}",content.1);
                    
                    kv = tokv(&mmap[content.0..=content.1]);
                    outkv.push(kv);
                }
            }

            let delzero:Vec<Vec<(i32, (usize, usize))>> = outkv.into_iter().map(
                |inner_vector| {
                    inner_vector
                        .into_iter()
                        .filter(|(i, _)| *i != 0)
                        .collect()
                }).collect();
                let mut outputlines:Vec<String> = Vec::new();
            outputlines = outputline(delzero, &title);
            
            let mut total_result = shared_result.lock().unwrap();
            total_result.push(outputlines);
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

    let mut file = File::create(cli.output).expect("unable to create file");
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
        "@ColNum 8".to_string(),
        "@ChrID 1".to_string(),
        "@ChrPos 2".to_string(),
        format!("@AC_{} 3", input_file_name),
        format!("@AG_{} 4", input_file_name),
        format!("@AT_{} 5", input_file_name),
        format!("@CA_{} 6", input_file_name),
        format!("@CG_{} 7", input_file_name),
        format!("@CT_{} 8", input_file_name),
        format!("@GA_{} 9", input_file_name),
        format!("@GC_{} 10", input_file_name),
        format!("@GT_{} 11", input_file_name),
        format!("@TA_{} 12", input_file_name),
        format!("@TC_{} 13", input_file_name),
        format!("@TG_{} 14", input_file_name),
        format!("@ins_{} 15", input_file_name),
        format!("@del_{} 16", input_file_name),
        format!("@ChrID_index {}",formatted_titles),
        format!("@ChrID_Strand {}", strand_info),
    ];
    for header in headers {
        writeln!(file, "{}", header).expect("Unable to write header");
    }
    for result in final_result.iter() {
        for line in result {
            writeln!(file,"{}",line).expect("unable to write value")
        }
    }
    let duration = start.elapsed();
    println!("Total run time: {:?}",duration);
}



fn tokv(slice: &[u8])->Vec<(i32,(usize,usize))>{
    let stringslice = String::from_utf8_lossy(slice).to_string();
    let shead:Vec<&str> = stringslice.split("\t").collect();
    let stail:String = shead[1].trim().to_string();
    let ssplit:Vec<&str> = stail.split(',').collect();
    let mut kvlist:Vec<(i32,(usize,usize))> = vec![];
    let mut sign:usize = 0;
    let mut part:Vec<&str>;
    let mut a:i32;
    let mut b:usize;

    for each in ssplit {
        if each.contains('x') {
            part = each.split('x').collect();
            a = part[0].parse::<i32>().expect("fail to i32");
            b = part[1].parse::<usize>().expect("fail to usize") - 1 + sign;
            kvlist.push((a,(sign, b)));
            sign = b + 1;
        } else {
            kvlist.push((each.parse::<i32>().expect("one word fail to i32"),(sign,sign)));
            sign = sign + 1;
        }
    }
    kvlist
}

fn outputline(input:Vec<Vec<(i32, (usize, usize))>>,title:&str) -> Vec<String>{
    let mut intervals: Vec<(usize, usize)> = Vec::new();
    for inner_vector in &input {
        for &(_, (start, end)) in inner_vector {
            intervals.push((start, end));
        }
    }
    intervals.sort_by(|a, b| a.0.cmp(&b.0)); 
    let mut merged_intervals: Vec<(usize, usize)> = Vec::new();
    for interval in intervals {
        if let Some(last) = merged_intervals.last_mut() {
            if interval.0 <= last.1 {
                last.1 = last.1.max(interval.1);
            } else {
                merged_intervals.push(interval);
            }
        } else {
            merged_intervals.push(interval);
        }
    }
    let mut output_vector: Vec<(usize, Vec<i32>)> = Vec::new();

    for &(start, end) in &merged_intervals {
        for i in start..=end {
            let mut row: Vec<i32> = vec![0; input.len()];

            for (col_idx, inner_vector) in input.iter().enumerate() {
                for &(value, (start, end)) in inner_vector {
                    if i >= start && i <= end {
                        row[col_idx] = value;
                    }
                }
            }
            output_vector.push((i, row));
        }
    }
    let mut formatted_output: Vec<String> = Vec::new();
    for (index, row) in output_vector.iter() {
        let row_str = row.iter().map(|v| v.to_string()).collect::<Vec<_>>().join(", ");
        let formatted_row = format!("{}, {}, {}", title.trim(), index + 1, row_str);
        formatted_output.push(formatted_row);
    }
    formatted_output
}

