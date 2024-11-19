use rayon::prelude::*;
use clap::Parser;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;
use std::sync::{Arc, Mutex};
use std::time::Instant;

#[derive(Debug, Clone)]
struct GeneEntry {
    chr_id: String,
    strand: char,
    position: u32,
    base1: Option<char>,
    base3: Option<char>,
    rt_1: Option<i32>,
    bd_1: Option<i32>,
    mutations: Vec<i32>,
    rt_3: Option<i32>,
    bd_3: Option<i32>,         
}

const BAR_LAB: &str = "-\\|/";

fn parse_file1(file_path: &str) -> io::Result<HashMap<(String, char, u32), GeneEntry>> {
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);

    let map: HashMap<(String, char, u32), GeneEntry> = reader
        .lines()
        .par_bridge() 
        .filter_map(|line| {
            let line = line.ok()?;
            if line.starts_with('@') {
                return None; 
            }
            let parts: Vec<&str> = line.split(',').map(|s| s.trim()).collect();
            if parts.len() < 5 {
                return None;
            }

            let chr_id = parts[0].to_string();
            let position: u32 = parts[1].parse().ok()?;
            let base1 = parts[2].chars().next();
            let rt_1: i32 = parts[3].parse().ok()?;
            let bd_1: i32 = parts[4].parse().ok()?;

            let entry = GeneEntry {
                chr_id: chr_id.clone(),
                strand: '+',
                position,
                base1,
                base3: None,
                rt_1: Some(rt_1),
                bd_1: Some(bd_1),
                mutations: vec![0; 14],
                rt_3: None,
                bd_3: None,
            };

            Some(((chr_id, '+', position), entry))
        })
        .collect();

    Ok(map)
}


fn parse_file2(file_path: &str) -> io::Result<HashMap<(String, char, u32), Vec<i32>>> {
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);

    let map: HashMap<(String, char, u32), Vec<i32>> = reader
        .lines()
        .par_bridge()
        .filter_map(|line| {
            let line = line.ok()?;
            if line.starts_with('@') {
                return None;
            }
            let parts: Vec<&str> = line.split(',').map(|s| s.trim()).collect();
            if parts.len() < 16 {
                return None;
            }

            let chr_id = parts[0].to_string();
            let position: u32 = parts[1].parse().ok()?;
            let mutations: Vec<i32> = parts[2..16]
                .iter()
                .map(|&s| s.parse().unwrap_or(0))
                .collect();

            Some(((chr_id, '+', position), mutations))
        })
        .collect();

    Ok(map)
}

// 解析文件3
fn parse_file3(file_path: &str) -> io::Result<HashMap<(String, char, u32), (Option<char>, i32, i32)>> {
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);

    let mut chr_id_index = vec![];
    let mut strand_index = vec![];

    let lines: Vec<_> = reader.lines().collect::<Result<_, _>>()?;

    for line in &lines {
        if line.starts_with("@ChrID_Index") {
            chr_id_index = line.split_whitespace().skip(1).map(String::from).collect();
            continue;
        }
        if line.starts_with("@ChrID_Strand") {
            strand_index = line.split_whitespace().skip(1).map(|s| s.chars().next().unwrap()).collect();
            continue;
        }
    }

    let map: HashMap<(String, char, u32), (Option<char>, i32, i32)> = lines
        .into_par_iter() 
        .filter_map(|line| {
            if line.starts_with('@') {
                return None; 
            }

            let parts: Vec<&str> = line.split(',').map(|s| s.trim()).collect();
            if parts.len() < 5 {
                return None; // 跳过不完整行
            }

            let index: usize = parts[0].parse::<usize>().ok()? - 1; // 减1以匹配索引
            let chr_id = chr_id_index[index].clone(); // 从索引中获取染色体ID
            let strand = strand_index[index]; // 从索引中获取链信息
            let position: u32 = parts[1].parse().ok()?;
            let base3 = parts[2].chars().next(); // 文件3的碱基信息
            let rt_3: i32 = parts[3].parse().ok()?;
            let bd_3: i32 = parts[4].parse().ok()?;

            Some(((chr_id, strand, position), (base3, rt_3, bd_3)))
        })
        .collect();

    Ok(map)
}

// 合并数据
fn merge_data(
    file1_data: HashMap<(String, char, u32), GeneEntry>,
    file2_data: HashMap<(String, char, u32), Vec<i32>>,
    file3_data: HashMap<(String, char, u32), (Option<char>, i32, i32)>
    ) -> Vec<GeneEntry> {
    let merged_entries = Arc::new(Mutex::new(HashMap::new())); // 使用Arc和Mutex包装HashMap
    let clen = &file1_data.len();
    let count = Arc::new(Mutex::new(0));
    // 使用par_iter并行处理
    file1_data.into_par_iter().for_each(|((chr_id, strand, position), entry)| {
        let counter = Arc::clone(&count);
        let mut merged_entries = merged_entries.lock().unwrap();
        merged_entries.insert((chr_id.clone(), strand, position), entry);
        let mut num = counter.lock().unwrap();
        *num += 1;
        print!(
            "\rmergefile1{}{}{}%",
            "#".repeat(((*num * 100) / clen)/ 2),
            BAR_LAB.chars().nth(((*num * 100) / clen) % 4).unwrap(),
            ((*num * 100) / clen)
        );
    });
    let count = Arc::new(Mutex::new(0));
    let clen = &file2_data.len();
    file2_data.into_par_iter().for_each(|((chr_id, strand, position), mutations)| {
        let counter = Arc::clone(&count);
        let mut merged_entries = merged_entries.lock().unwrap();
        let entry = merged_entries.entry((chr_id.clone(), strand, position)).or_insert(GeneEntry {
            chr_id: chr_id.clone(),
            strand,
            position,
            base1: None,
            base3: None,
            rt_1: None,
            bd_1: None,
            mutations: vec![0; 14],
            rt_3: None,
            bd_3: None,
        });

        entry.mutations = mutations;
        let mut num = counter.lock().unwrap();
        *num += 1;
        print!(
            "\rmergefile2{}{}{}%",
            "#".repeat(((*num * 100) / clen)/ 2),
            BAR_LAB.chars().nth(((*num * 100) / clen) % 4).unwrap(),
            ((*num * 100) / clen)
        );
    });

    let count = Arc::new(Mutex::new(0));
    let clen = &file3_data.len();
    file3_data.into_par_iter().for_each(|((chr_id, strand, position), (base3, rt_3, bd_3))| {
        let counter = Arc::clone(&count);
        let mut merged_entries = merged_entries.lock().unwrap();
        let entry = merged_entries.entry((chr_id.clone(), strand, position)).or_insert(GeneEntry {
            chr_id: chr_id.clone(),
            strand,
            position,
            base1: None,
            base3,
            rt_1: None,
            bd_1: None,
            mutations: vec![0; 14],
            rt_3: None,
            bd_3: None,
        });
        
        entry.base3 = base3.or(entry.base3);
        entry.rt_3 = Some(rt_3);
        entry.bd_3 = Some(bd_3);
        let mut num = counter.lock().unwrap();
        *num += 1;
        print!(
            "\rmergefile3{}{}{}%",
            "#".repeat(((*num * 100) / clen)/ 2),
            BAR_LAB.chars().nth(((*num * 100) / clen) % 4).unwrap(),
            ((*num * 100) / clen)
        );
    });

    let mut result: Vec<GeneEntry> = merged_entries.lock().unwrap().clone().into_iter().map(|(_, v)| v).collect();
    
    result.par_sort_by(|a, b| {
        a.chr_id.cmp(&b.chr_id)
            .then_with(|| a.strand.cmp(&b.strand))
            .then_with(|| a.position.cmp(&b.position))
    });
    result
}

// 写入到新的CSV文件
fn write_to_csv(file_path: &str, data: Vec<GeneEntry>) -> io::Result<()> {
    let path = Path::new(file_path);
    let mut file = File::create(&path)?;

    // 写入表头
    writeln!(file, "ChrID,Strand,Position,Base1,RT1,BD1,AC,AG,AT,CA,CG,CT,GA,GC,GT,TA,TC,TG,Ins,Del,Base3,RT3,BD3")?;

    for entry in data {
        // 如果base1或base3为空，用'N'代替
        let base1_str = entry.base1.map_or("N".to_string(), |b| b.to_string());
        let base3_str = entry.base3.map_or("N".to_string(), |b| b.to_string());


        // 将每个条目写入CSV文件
        writeln!(
            file,
            "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
            entry.chr_id,
            entry.strand,
            entry.position,
            base1_str,
            entry.rt_1.unwrap_or(0),
            entry.bd_1.unwrap_or(0),
            entry.mutations[0], entry.mutations[1], entry.mutations[2], entry.mutations[3],
            entry.mutations[4], entry.mutations[5], entry.mutations[6], entry.mutations[7],
            entry.mutations[8], entry.mutations[9], entry.mutations[10], entry.mutations[11],
            entry.mutations[12], entry.mutations[13],
            base3_str, // 文件3的碱基信息
            entry.rt_3.unwrap_or(0),
            entry.bd_3.unwrap_or(0)
        )?;
    }

    Ok(())
}

#[derive(Parser)]
#[command(name="merge", author="hyf", version="1.0", about="merge three files", long_about = None)]
struct Cli {
    #[arg(short,long)]
    csv: String,
    #[arg(short,long)]
    pipe: String,
    #[arg(short,long)]
    txt: String,
    #[arg(short,long)]
    output: String,
}


fn main() -> io::Result<()> {
    let now = Instant::now();
    let cli = Cli::parse();
    let file_path1 = &cli.csv;
    let file_path2 = &cli.pipe;
    let file_path3 = &cli.txt;

    // 使用 rayon::join 来并行运行三个任务
    let (result1, result2, result3) = {
        let mut result1 = None;
        let mut result2 = None;
        let mut result3 = None;

        rayon::scope(|s| {
            s.spawn(|_| {
                result1 = Some(parse_file1(file_path1));
            });

            s.spawn(|_| {
                result2 = Some(parse_file2(file_path2));
            });

            s.spawn(|_| {
                result3 = Some(parse_file3(file_path3));
            });
        });

        (result1, result2, result3)
    };

    // 确保 result 变量在 scope 结束后被 unwrap
    let file1_data = result1.unwrap()?;
    let file2_data = result2.unwrap()?;
    let file3_data = result3.unwrap()?;



    println!("three file read over");
    let output_file = &cli.output;
    println!("merge data");
    let merged_data = merge_data(file1_data, file2_data, file3_data);
    println!("write to csv");
    write_to_csv(output_file, merged_data)?;

    let end = now.elapsed().as_secs();
    
    println!("Total runtime {:?} s",end);
    Ok(())
}
