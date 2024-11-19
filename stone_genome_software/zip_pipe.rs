use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, Result, Write};

fn main() -> Result<()> {
    // Read command line arguments
    let args: Vec<String> = env::args().collect();
    if args.len() != 3 {
        eprintln!("Usage: {} <input_file> <output_file>", args[0]);
        std::process::exit(1);
    }

    let input_file = &args[1];
    let output_file = &args[2];

    let file = File::open(input_file)?;
    let reader = BufReader::new(file);

    let mut chrid_strand_list = Vec::new();
    let mut file_header = Vec::new();
    let mut optimized_lines = Vec::new();
    let mut strand_info = Vec::new(); // New vector to store strand information

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('@') {
            file_header.push(line.clone()); // Clone the line for the header
            if line.starts_with("@ChrID_Index") {
                file_header.push("@ChrID_Strand".to_string()); // Add @ChrID_Strand after @ChrID_Index
            }
        } else {
            let parts: Vec<String> = line.split(',').map(|s| s.to_string()).collect();

            let chrid = &parts[0];
            let strand = &parts[1];
            let chrid_strand = format!("{}({})", chrid, strand);

            if !chrid_strand_list.contains(&chrid_strand) {
                chrid_strand_list.push(chrid_strand.clone());
                strand_info.push(strand.to_string()); // Store strand information
            }

            let chrid_index = chrid_strand_list.iter().position(|r| r == &chrid_strand).unwrap() + 1;
            let chrid_index_str = chrid_index.to_string();

            let mut new_parts = parts.clone();
            new_parts[0] = chrid_index_str;
            new_parts.remove(1); // Remove the strand information
            optimized_lines.push(new_parts.join(","));
        }
    }

    let mut optimized_file = File::create(output_file)?;

    for line in &file_header {
        writeln!(optimized_file, "{}", line)?;
    }

    // Write the index section
    writeln!(optimized_file, "@ChrID_Index\t{}", chrid_strand_list.iter().map(|s| s.split('(').next().unwrap()).collect::<Vec<&str>>().join("\t"))?;

    // Write the strand information section
    writeln!(optimized_file, "@ChrID_Strand\t{}", strand_info.join("\t"))?;

    // Write the optimized data lines
    for line in optimized_lines {
        writeln!(optimized_file, "{}", line)?;
    }

    Ok(())
}

