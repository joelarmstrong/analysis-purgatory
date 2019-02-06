use std::collections::{HashSet, HashMap};
use std::env;
use std::fmt;
use std::fs::File;
use std::io;
use std::io::{Read, BufReader, Write};
use std::process::exit;
use failure::{Fallible, bail, ensure, ResultExt};
use csv;

fn main() {
    let argv: Vec<String> = env::args().collect();
    let subcommand = argv.get(1).expect("Must have 2 arguments (find/fix and c2h file to check)");
    let path = argv.get(2).expect("Must have 2 arguments (find/fix and c2h file to check)");
    let reader = BufReader::new(File::open(path).expect("Can't open input file"));
    let chrom_sizes_path = argv.get(3).expect("Must provide chrom.sizes file as 3rd argument");
    let chrom_sizes = read_chrom_sizes(chrom_sizes_path).expect("couldn't read chrom.sizes");
    match subcommand.as_ref() {
        "find" => find(reader, &chrom_sizes, path),
        "fix" => {
            fix(reader, &chrom_sizes, &mut io::stdout().lock()).expect("fix failed")
        },
        _ => panic!("Only find and fix are valid subcommnands"),
    }
}

fn find(input: impl Read, chrom_sizes: &HashMap<String, u64>, path: &str) {
    match find_buggy_lines(input, chrom_sizes) {
        Ok(bugs) => {
            if !bugs.is_empty() {
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

fn read_chrom_sizes(path: &str) -> Fallible<HashMap<String, u64>> {
    let mut chrom_sizes = HashMap::new();
    let rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_path(path)?;
    for result in rdr.into_records() {
        let record = result?;
        chrom_sizes.insert(record[0].into(), record[1].parse()?);
    }
    Ok(chrom_sizes)
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

/// Fix a region which ends in a buggy region. Return the names of the
/// skipped bugged segments.
fn fix_and_output(buf: &mut Vec<C2HLine>, start_after_bad_region: u64, output: &mut impl Write) -> Fallible<Vec<u64>> {
    let mut skip: Vec<bool> = vec![false; buf.len()];
    // Iterate backwards through the problematic region. Bad lines
    // will have a length field containing the true start position of
    // the previous line, so they are easily detectable this way.
    let mut true_start = start_after_bad_region;
    for i in (0..buf.len()).rev() {
        match buf[i].data {
            TopSegment { start, length, .. } | BottomSegment { start, length, .. } => {
                if length == true_start && start != 0 {
                    // Bad line.
                    skip[i] = true;
                } else {
                    // Good line.
                    true_start -= length;
                }
            }
            _ => {},
        }
    }

    // Now that we know which lines are bad, we can go through in the
    // forward direction and only output good lines, recalculating the
    // true start position ourselves.
    for i in (0..buf.len()).filter(|&i| !skip[i]) {
        match &mut buf[i].data {
            TopSegment { start, length, .. } | BottomSegment { start, length, .. } => {
                *start = true_start;
                true_start += *length;
            }
            _ => {},
        }
        write!(output, "{}", buf[i].data)?;
    }
    let removed = buf.iter()
        .enumerate()
        .filter(|(i, _)| skip[*i])
        .map(|(_, l)| match l.data {
            TopSegment { name, .. } => {
                name.unwrap_or(0)
            },
            BottomSegment { name, .. } => {
                name
            },
            _ => unreachable!(),
        })
        .collect();
    buf.clear();
    Ok(removed)
}

fn fix(input: impl Read, chrom_sizes: &HashMap<String, u64>, output: &mut impl Write) -> Fallible<()> {
    let mut buf: Vec<C2HLine> = vec![];
    let mut prev_end = 0;
    let mut removed_segments = HashSet::new();
    let mut cur_seq_length = None;
    for result in get_iter(input)? {
        let c2h_line = result?;
        match &c2h_line.data {
            BottomSegment { start, length, .. } | TopSegment { start, length, .. } => {
                if *start < prev_end {
                    // End of problem region. Go back through and fix what needs to be fixed.
                    println!("problem region found at {}", *start);
                    removed_segments.extend(fix_and_output(&mut buf, *start, output)?);
                }
                prev_end = start + length;
            },
            NewSequence { header, .. } => {
                if let Some(&seq_length) = cur_seq_length {
                    if seq_length != prev_end {
                        removed_segments.extend(fix_and_output(&mut buf, seq_length, output)?);
                    }
                }
                for item in buf.drain(..) {
                    write!(output, "{}", item.data)?;
                }
                prev_end = 0;
                cur_seq_length = chrom_sizes.get(&header[1..header.len() - 1]);
            }
        }
        if let TopSegment { name: Some(name), .. } = c2h_line.data {
            // Ensure this isn't referencing a removed bottom segment.
            ensure!(!removed_segments.contains(&name),
                    "Top segment referenced removed bottom segment {}", name);
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

fn find_buggy_lines(input: impl Read, chrom_sizes: &HashMap<String, u64>) -> Fallible<Vec<u64>> {
    let mut ret = vec![];
    let mut cur_pos = 0;
    let mut cur_seq_length: Option<&u64> = None;
    for result in get_iter(input)? {
        let c2h = result?;
        match c2h.data {
            TopSegment { start, length, .. } | BottomSegment { start, length, .. } => { 
                if cur_pos != start {
                    ret.push(c2h.line_num);
                }
                cur_pos = start + length;
            },
            NewSequence { header, .. } => {
                if let Some(cur_seq_length) = cur_seq_length {
                    if cur_pos != *cur_seq_length {
                        ret.push(c2h.line_num);
                    }
                }
                cur_seq_length = chrom_sizes.get(&header);
                cur_pos = 0;
            },
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
        fix(input.as_bytes(), &HashMap::new(), &mut output).unwrap();
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
        fix(input.as_bytes(), &HashMap::new(), &mut output).unwrap();
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
        fix(input.as_bytes(), &HashMap::new(), &mut output).unwrap();
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

        // OK, this has to be the last way that things can get messed
        // up
        let input = "s	'birdAnc248'	'birdAnc248refChr2611'	1
a	2633057669639824667	0	6
a	8856750879689135289	6	4
a	2633057669639824664	10	5
a	7770820419539392828	15	18
a	5176324821708937241	33	7
a	5176324821708937239	40	40
a	4183562578850455601	80	1
a	4183562578850455599	81	41
a	4183562578850455598	122	7
a	4183562578850455596	129	48
a	5176324821708937238	48	7
a	7571113923563180385	55	54
a	5176324821708937235	109	6
a	5176324821708937232	115	7
a	7571113923563180382	122	40
a	5176324821708937229	162	8
a	5176324821708937227	170	170
a	4183562578850455595	340	67
a	4183562578850455593	407	237
a	4183562578850455592	644	25
a	4183562578850455590	669	262
a	4183562578850455589	931	339
a	4183562578850455587	1270	601
a	4183562578850455586	1871	26
a	4183562578850455584	1897	627
a	4183562578850455583	2524	19
a	4183562578850455581	2543	646
a	4183562578850455580	3189	31
a	4183562578850455578	3220	677
a	4183562578850455577	3897	98
a	4183562578850455575	3995	775
a	4183562578850455574	4770	55
a	4183562578850455572	4825	830
a	4183562578850455571	5655	1
a	4183562578850455569	5656	831
a	5176324821708937226	831	12
a	7571113923563180379	843	351";
        let mut output = vec![];
        fix(input.as_bytes(), &HashMap::new(), &mut output).unwrap();
        assert_eq!(from_utf8(&output).unwrap(), "s	'birdAnc248'	'birdAnc248refChr2611'	1
a	2633057669639824667	0	6
a	8856750879689135289	6	4
a	2633057669639824664	10	5
a	7770820419539392828	15	18
a	5176324821708937241	33	7
a	4183562578850455601	40	1
a	4183562578850455598	41	7
a	5176324821708937238	48	7
a	7571113923563180385	55	54
a	5176324821708937235	109	6
a	5176324821708937232	115	7
a	7571113923563180382	122	40
a	5176324821708937229	162	8
a	4183562578850455595	170	67
a	4183562578850455592	237	25
a	4183562578850455589	262	339
a	4183562578850455586	601	26
a	4183562578850455583	627	19
a	4183562578850455580	646	31
a	4183562578850455577	677	98
a	4183562578850455574	775	55
a	4183562578850455571	830	1
a	5176324821708937226	831	12
a	7571113923563180379	843	351
");

        let input = "s	'birdAnc64'	'birdAnc64refChr210'	1
a	7663859928389344295	0	142
a	3460734838657496896	142	24
a	3460734838657496893	166	6
a	3460734838657496890	172	48
a	3460734838657496887	220	160
a	3460734838657496884	380	126
a	3460734838657496881	506	33
a	3460734838657496878	539	98
a	3460734838657496875	637	41
a	3460734838657496872	678	11
a	7663859928389344298	689	9
a	1041598151317810649	698	5
a	1041598151317810652	703	5
a	1930073915304958345	708	18
a	1041598151317810646	726	6
a	1041598151317810644	732	732
a	4197777065174346908	1464	18
a	4197777065174346906	1482	750
a	1041598151317810643	750	7
a	1930073915304958342	757	90
a	1041598151317810640	847	4
a	1041598151317810638	851	851
a	4197777065174346905	1702	2
a	4197777065174346903	1704	853
a	1041598151317810637	853	2
a	1041598151317810634	855	16
a	1930073915304958339	871	143
a	1041598151317810631	1014	23
a	1041598151317810628	1037	32
a	1041598151317810626	1069	1069
a	4197777065174346902	2138	10
a	4197777065174346900	2148	1079
a	1041598151317810625	1079	2
a	1041598151317810623	1081	1081
a	4197777065174346899	2162	1
a	4197777065174346897	2163	1082
a	1041598151317810622	1082	39
a	1041598151317810620	1121	1121
a	4197777065174346896	2242	2
a	4197777065174346894	2244	1123
a	1041598151317810619	1123	13
a	1041598151317810617	1136	1136
a	4197777065174346893	2272	1
a	4197777065174346891	2273	1137
a	1041598151317810616	1137	13
a	1041598151317810614	1150	1150
a	4197777065174346890	2300	1
a	4197777065174346888	2301	1151
a	1041598151317810613	1151	109
a	1041598151317810610	1260	17
a	1930073915304958336	1277	40
a	1041598151317810607	1317	26
a	1930073915304958333	1343	134
a	1041598151317810604	1477	14
a	1930073915304958330	1491	24
a	1930073915304958328	1515	1515
a	4197777065174346873	3030	25
a	4197777065174346871	3055	1540
a	4197777065174346876	4595	4
a	5245849140956493994	4599	38
a	5245849140956493991	4637	49
a	4197777065174346879	6143	13
a	4197777065174346877	6156	1644
a	4197777065174346882	7800	91
a	4197777065174346880	7891	1735
s	'blah'	'foo'	1";
        let mut output = vec![];
        let mut map = HashMap::new();
        map.insert("birdAnc64refChr210".to_owned(), 1735);
        fix(input.as_bytes(), &map, &mut output).unwrap();
        assert_eq!(from_utf8(&output).unwrap(), "");
    }

    /// Ensure a lost segment that is used later on causes a crash.
    #[test]
    fn test_fix_paranoia() {
        let input = "s	'parent'	'parentChr'	1
a	0	0	366
a	4799148352916656454	366	14
a	4799148352916656452	380	380
a	4183562578850463929	760	1
a	4183562578850463927	761	381
a	4799148352916656457	381	2
s	'child'	'childChr'	0
a	0	1234	8675309	1
a	1234	380	4799148352916656452	1
a	1614	1	4183562578850463929	1
";
        let mut output = vec![];
        let result = fix(input.as_bytes(), &HashMap::new(), &mut output);
        assert!(result.is_err());
    }
}
