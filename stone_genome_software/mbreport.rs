use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufReader, BufRead, Write, BufWriter};
use clap::Parser;

#[derive(Parser)]
#[command(name="mbreport", author="hyf", version="1.0", about="statistic", long_about = None)]
struct Cli {
    #[arg(short,long)]
    input: String,
    #[arg(short,long)]
    depth: u32,
    #[arg(short,long)]
    output: String,
}

fn main() -> io::Result<()> {
    let cli = Cli::parse();
    let file = File::open(cli.input)?; 
    let reader = BufReader::new(file); 

    let mut data_map: HashMap<(String, String, String), (u32, u32, u32, u32)> = HashMap::new();
    let lines = reader.lines(); 

    for line in lines {
        let line = line?; 
        if line.trim().is_empty() {
            continue; 
        }
        let fields: Vec<&str> = line.split(',').collect(); 
        let chrid = fields[0].to_string();
        let geneid = fields[1].to_string();
        let position = fields[3].to_string();

        let rf_mutation_count: u32 = fields[6].parse().unwrap_or(0);
        let rf_mutation_depth: u32 = fields[7].parse().unwrap_or(0);
        let pipe_truncation_count: u32 = fields[23].parse().unwrap_or(0);
        let pipe_truncation_bd: u32 = fields[24].parse().unwrap_or(0);


        if rf_mutation_depth > cli.depth && pipe_truncation_bd > cli.depth && rf_mutation_count*4 < rf_mutation_depth && pipe_truncation_count*4 < pipe_truncation_bd{
            data_map.insert((chrid, geneid, position), (rf_mutation_count, rf_mutation_depth, pipe_truncation_count, pipe_truncation_bd));
        }
}
    let infomap = calculate_gene_statistics(data_map);

    let file = File::create(cli.output)?;
    let mut writer = BufWriter::new(file);


    writeln!(writer, "GeneID,Avg_RF_Count,Avg_RF_Depth,RF_Ratio,Avg_Pipe_Count,Avg_Pipe_BD,Pipe_Ratio")?;


    for (geneid, (avg_rf_count, avg_rf_depth, rf_ratio, avg_pipe_count, avg_pipe_bd, pipe_ratio)) in &infomap {
        writeln!(
            writer,
            "{},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6}",
            geneid, avg_rf_count, avg_rf_depth, rf_ratio, avg_pipe_count, avg_pipe_bd, pipe_ratio
        )?;
    }

  
    writer.flush()?;

    Ok(())
}

fn calculate_gene_statistics(
    data_map: HashMap<(String, String, String), (u32, u32, u32, u32)>
) -> HashMap<String, (f64, f64, f64, f64, f64, f64)> {
    let mut gene_stats_map: HashMap<String, (f64, f64, f64, f64, f64, f64)> = HashMap::new();
    let mut gene_aggregates: HashMap<String, (u32, u32, u32, u32, usize)> = HashMap::new();


    for ((_, geneid, _), (rf_count, rf_depth, pipe_count, pipe_bd)) in &data_map {
        let entry = gene_aggregates.entry(geneid.clone()).or_insert((0, 0, 0, 0, 0));
        entry.0 += *rf_count;
        entry.1 += *rf_depth;
        entry.2 += *pipe_count;
        entry.3 += *pipe_bd;
        entry.4 += 1; //
    }


    for (geneid, (total_rf_count, total_rf_depth, total_pipe_count, total_pipe_bd, count)) in gene_aggregates {
        let avg_rf_count = total_rf_count as f64 / count as f64;
        let avg_rf_depth = total_rf_depth as f64 / count as f64;
        let rf_ratio = avg_rf_count / avg_rf_depth;
        let avg_pipe_count = total_pipe_count as f64 / count as f64;
        let avg_pipe_bd = total_pipe_bd as f64 / count as f64;
        let pipe_ratio = avg_pipe_count / avg_pipe_bd;

        gene_stats_map.insert(geneid, (avg_rf_count, avg_rf_depth, rf_ratio, avg_pipe_count, avg_pipe_bd, pipe_ratio));
    }

    gene_stats_map
}