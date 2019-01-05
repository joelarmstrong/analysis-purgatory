use std::fmt;
use std::io;
use std::io::{Read, BufReader, Write};
use std::fs::File;
use std::env;
use std::process::exit;
use failure::{Fallible, bail, ResultExt};
use csv;

fn main() {
    let argv: Vec<String> = env::args().collect();
    let subcommand = argv.get(1).expect("Must have 2 arguments (find/fix and c2h file to check)");
    let path = argv.get(2).expect("Must have 2 arguments (find/fix and c2h file to check)");
    let reader = BufReader::new(File::open(path).expect("Can't open input file"));
    match subcommand.as_ref() {
        "find" => find(reader, path),
        "fix" => fix(reader, &mut io::stdout().lock()).expect("fix failed"),
        _ => panic!("Only find and fix are valid subcommnands"),
    }
}

fn find(input: impl Read, path: &str) {
    match find_buggy_lines(input) {
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

enum C2HData {
    BottomSegment { name: u64, start: u64, length: u64 },
    TopSegment { name: Option<u64>, start: u64, length: u64, reversed: Option<u8> },
    NewSequence { event: String, header: String, is_bottom: u8 },
}

use self::C2HData::*;

struct C2HLine {
    line_num: u64,
    data: C2HData,
}

impl fmt::Display for C2HData {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match &self {
            BottomSegment { name, start, length } => {
                writeln!(f, "a\t{}\t{}\t{}", name, start, length)
            },
            TopSegment { name, start, length, reversed } => {
                match (name, reversed) {
                    (Some(n), Some(r)) => writeln!(f, "a\t{}\t{}\t{}\t{}", start, length, n, r),
                    (None, None) => writeln!(f, "a\t{}\t{}", start, length),
                    _ => panic!("Got name or reversed but not both in top segment"),
                }
            },
            NewSequence { event, header, is_bottom } => {
                writeln!(f, "s\t{}\t{}\t{}", event, header, is_bottom)
            }
        }
    }
}

fn fix_and_output(buf: &mut Vec<C2HLine>, start_after_bad_region: u64, output: &mut impl Write) -> Fallible<()> {
    println!("fixing region that ends in good start pos {}", start_after_bad_region);
    let mut in_region = false;
    let mut next_pos = 0;
    let mut prev_end = 0;
    for mut item in buf.drain(..) {
        match &mut item.data {
            TopSegment { start: other_start, length: other_length, .. } | BottomSegment { start: other_start, length: other_length, .. } => {
                if !in_region && *other_start + *other_length >= start_after_bad_region {
                    // Beginning of problem region.
                    in_region = true;
                }
                if in_region && next_pos != 0 {
                    *other_start = next_pos;
                    next_pos = 0;
                    prev_end = *other_start + *other_length;
                    write!(output, "{}", item.data)?;
                } else if in_region {
                    if *other_length == prev_end {
                        // Extraneous, skip
                        next_pos = *other_length;
                    } else {
                        *other_start = prev_end;
                        prev_end = *other_start + *other_length;
                        write!(output, "{}", item.data)?;
                    }
                } else {
                    prev_end = *other_start + *other_length;
                    write!(output, "{}", item.data)?;
                }
            }
            NewSequence { .. } => {
                write!(output, "{}", item.data)?;
                prev_end = 0;
            },
        }
    }
    Ok(())
}

fn fix(input: impl Read, output: &mut impl Write) -> Fallible<()> {
    let mut buf: Vec<C2HLine> = vec![];
    let mut prev_start = 0;
    for result in get_iter(input)? {
        let c2h_line = result?;
        match c2h_line.data {
            BottomSegment { start, .. } | TopSegment { start, .. } => {
                if start < prev_start {
                    // End of problem region. Go back through and fix what needs to be fixed.
                    fix_and_output(&mut buf, start, output)?;
                }
                prev_start = start;
            },
            NewSequence { .. } => {
                prev_start = 0;
                for item in buf.drain(..) {
                    write!(output, "{}", item.data)?;
                }
            }
        }
        buf.push(c2h_line);
    }
    for item in buf.into_iter() {
        write!(output, "{}", item.data)?;
    }
    Ok(())
}

fn get_iter(input: impl Read) -> Fallible<impl Iterator<Item = Fallible<C2HLine>>> {
    let rdr = csv::ReaderBuilder::new()
        .flexible(true)
        .has_headers(false)
        .delimiter(b'\t')
        .from_reader(input);
    Ok(rdr.into_records().map(|result| -> Fallible<C2HLine> {
        let record = result?;
        let line = record.position().map(|p| p.line()).unwrap_or(0);
        if &record[0] == "a" {
            if record.len() == 4 {
                let name = record[1].parse::<u64>().with_context(|_| format!("reading name at line {}, line", line))?;
                let start: u64 = record[2].parse::<u64>().with_context(|_| format!("reading position at line {}", line))?;
                let length: u64 = record[3].parse::<u64>().with_context(|_| format!("reading length at line {}", line))?;
                Ok(C2HLine { data: BottomSegment {name, start, length }, line_num: line })
            } else {
                let start: u64 = record[1].parse::<u64>().with_context(|_| format!("reading position at line {}", line))?;
                let length: u64 = record[2].parse::<u64>().with_context(|_| format!("reading length at line {}", line))?;
                let name: Option<u64> = record.get(3).and_then(|s| s.parse::<u64>().ok());
                let reversed: Option<u8> = record.get(4).and_then(|s| s.parse::<u8>().ok());
                Ok(C2HLine { data: TopSegment { start, length, name, reversed }, line_num: line })
            }
        } else if &record[0] == "s" {
            let event = record[1].to_owned();
            let header = record[2].to_owned();
            let is_bottom = record[3].parse::<u8>()?;
            Ok(C2HLine { data: NewSequence { event, header, is_bottom }, line_num: line })
        } else {
            bail!("failure to understand line {}", line);
        }
    }))
}

fn find_buggy_lines(input: impl Read) -> Fallible<Vec<u64>> {
    let mut ret = vec![];
    let mut cur_pos = 0;
    for result in get_iter(input)? {
        let c2h = result?;
        match c2h.data {
            TopSegment { start, length, .. } | BottomSegment { start, length, .. } => { 
                if cur_pos != start {
                    ret.push(c2h.line_num);
                }
                cur_pos = start + length;
            },
            NewSequence { .. } => cur_pos = 0,
        }
    }
    Ok(ret)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::str::from_utf8;
    use pretty_assertions::assert_eq;

    #[test]
    fn test_fix() {
        let input = "s	'foo'	'bar'	0
a	0	0	327
a	5572078638964173662	327	13
a	3360248271971874541	340	20
a	5572078638964173668	360	6
a	4799148352916656454	366	14
a	4799148352916656452	380	380
a	4183562578850463929	760	1
a	4183562578850463927	761	381
a	4799148352916656457	381	2
";
        let mut output = vec![];
        fix(input.as_bytes(), &mut output).unwrap();
        assert_eq!(from_utf8(&output).unwrap(), "s	'foo'	'bar'	0
a	0	0	327
a	5572078638964173662	327	13
a	3360248271971874541	340	20
a	5572078638964173668	360	6
a	4799148352916656454	366	14
a	4183562578850463929	380	1
a	4799148352916656457	381	2
");

        // Second test
        let input = "a	0	0	1020209
a	7043489079719062257	1020209	3
a	7043489079719062255	1020212	1020212
a	6533315684431011826	2040424	5
a	6533315684431011824	2040429	1020217
a	6533315684431011823	3060646	5
a	6533315684431011821	3060651	1020222
a	7043489079719062260	1020222	71
";
        let mut output = vec![];
        fix(input.as_bytes(), &mut output).unwrap();
        assert_eq!(from_utf8(&output).unwrap(), "a	0	0	1020209
a	7043489079719062257	1020209	3
a	6533315684431011826	1020212	5
a	6533315684431011823	1020217	5
a	7043489079719062260	1020222	71
");

        // Third test
        let input = "
a	2900458897514937132	0	43425
a	2900458897514937132	43425	8
a	2900458897514937135	43433	7
a	2900458897514937138	43440	1
a	2900458897514937136	43441	43441
a	3871266092189995034	86882	17
a	3871266092189995032	86899	43458
a	2900458897514937141	43458	7
a	2900458897514937139	43465	43465
a	3871266092189995031	86930	5
a	1582593056555652163	86935	3
a	3871266092189995028	130405	1
a	3871266092189995026	130406	43474
a	2900458897514937144	43474	6
a	2900458897514937142	43480	43480
a	3871266092189995025	86960	8
a	3871266092189995023	86968	43488
a	2900458897514937147	43488	23
a	2900458897514937150	43511	12
a	2900458897514937153	43523	13
a	2900458897514937156	43536	14
";
        let mut output = vec![];
        fix(input.as_bytes(), &mut output).unwrap();
        assert_eq!(from_utf8(&output).unwrap(), "a	2900458897514937132	0	43425
a	2900458897514937132	43425	8
a	2900458897514937135	43433	7
a	2900458897514937138	43440	1
a	3871266092189995034	43441	17
a	2900458897514937141	43458	7
a	3871266092189995031	43465	5
a	1582593056555652163	43470	3
a	3871266092189995028	43473	1
a	2900458897514937144	43474	6
a	3871266092189995025	43480	8
a	2900458897514937147	43488	23
a	2900458897514937150	43511	12
a	2900458897514937153	43523	13
a	2900458897514937156	43536	14
");
    }
}
