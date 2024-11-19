use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Seek, SeekFrom, Write};
use std::sync::{Arc, Mutex};
use std::thread;

fn process_line(line: &str) -> String {
    if let Some((_, rest)) = line.split_once('\t') {
        let parts: Vec<&str> = rest.split(',').collect();
        let mut compressed = Vec::new();
        let mut count = 1;

        for i in 1..parts.len() {
            if parts[i] == parts[i - 1] {
                count += 1;
            } else {
                if count > 1 {
                    compressed.push(format!("{}x{}", parts[i - 1], count));
                } else {
                    compressed.push(parts[i - 1].to_string());
                }
                count = 1;
            }
        }
        if count > 1 {
            compressed.push(format!("{}x{}", parts[parts.len() - 1], count));
        } else {
            compressed.push(parts[parts.len() - 1].to_string());
        }

        format!("{}\t{}", line.split_once('\t').unwrap().0, compressed.join(","))
    } else {
        line.to_string()
    }
}

fn main() -> io::Result<()> {
    let args: Vec<String> = env::args().collect();
    if args.len() != 4 {
        eprintln!("Usage: {} <input_file> <output_file> <num_threads>", args[0]);
        std::process::exit(1);
    }

    let input_path = &args[1];
    let output_path = &args[2];
    let num_threads: usize = args[3].parse().unwrap_or_else(|_| {
        eprintln!("Invalid number of threads");
        std::process::exit(1);
    });

    let input_file = File::open(input_path)?;
    let output_file = File::create(output_path)?;
    let writer = Arc::new(Mutex::new(BufWriter::new(output_file)));

    let file_size = input_file.metadata()?.len();
    let chunk_size = file_size / num_threads as u64;

    let mut handles = Vec::new();
    let result = Arc::new(Mutex::new(Vec::new()));

    for i in 0..num_threads {
        let input_path = input_path.clone();
        let result = Arc::clone(&result);

        let handle = thread::spawn(move || -> io::Result<()> {
            let mut input_file = File::open(input_path)?;
            input_file.seek(SeekFrom::Start(i as u64 * chunk_size))?;
            let mut reader = BufReader::new(input_file);

            let mut local_lines = Vec::new();
            let mut processed_size = 0;

            if i != 0 {
                // 如果不是第一个块，跳过第一个完整的行，以避免跨块重复
                let mut first_line = String::new();
                reader.read_line(&mut first_line)?;
                processed_size += first_line.len() as u64;
            }

            for line in reader.lines() {
                let line = line?;
                if line.trim().is_empty() {
                    continue;
                }
                local_lines.push((i as u64 * chunk_size + processed_size, process_line(&line)));
                processed_size += line.len() as u64;

                if processed_size >= chunk_size {
                    break;
                }
            }

            let mut result = result.lock().unwrap();
            result.extend(local_lines);

            Ok(())
        });
        handles.push(handle);
    }

    for handle in handles {
        handle.join().unwrap()?;
    }

    let mut result = result.lock().unwrap();
    result.sort_by_key(|&(idx, _)| idx);

    let mut writer = writer.lock().unwrap();
    for (_, line) in result.iter() {
        writeln!(writer, "{}", line)?;
    }

    Ok(())
}

