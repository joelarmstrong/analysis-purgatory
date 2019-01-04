use std::env;
use std::process::exit;
use failure::{Fallible, bail, ResultExt};
use csv;

fn main() {
    let argv: Vec<String> = env::args().collect();
    let path = argv.get(1).expect("Must have 1 argument (c2h file to check)");
    match run(path) {
        Ok(bugs) => {
            if bugs.len() > 0 {
                eprintln!("Found {} bugs in file {}, lines:", bugs.len(), path);
                for bug in bugs {
                    eprintln!("{}", bug);
                }
                exit(1);
            }
        }
        Err(e) => {
            panic!("Got error when parsing: {}", e);
        }
    }
}

fn run(path: &str) -> Fallible<Vec<u64>> {
    let mut ret = vec![];
    let mut rdr = csv::ReaderBuilder::new()
        .flexible(true)
        .has_headers(false)
        .delimiter(b'\t')
        .from_path(path)?;
    let mut cur_pos = 0;
    for result in rdr.records() {
        let record = result?;
        let line = record.position().map(|p| p.line()).unwrap_or(0);
        if &record[0] == "a" {
            if record.len() == 4 {
                let pos: u64 = record[2].parse::<u64>().with_context(|_| format!("reading position at line {}", line))?;
                let length: u64 = record[3].parse::<u64>().with_context(|_| format!("reading position at line {}", line))?;
                if cur_pos != pos {
                    ret.push(line);
                }
                cur_pos = pos + length;
            } else {
                let pos: u64 = record[1].parse::<u64>().with_context(|_| format!("reading position at line {}", line))?;
                let length: u64 = record[2].parse::<u64>().with_context(|_| format!("reading position at line {}", line))?;
                if cur_pos != pos {
                    ret.push(line);
                }
                cur_pos = pos + length;
            }
        } else if &record[0] == "s" {
            cur_pos = 0;
        } else {
            bail!("failure to understand line {}", line);
        }
    }
    Ok(ret)
}
