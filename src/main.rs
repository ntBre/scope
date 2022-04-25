use std::str::FromStr;
use std::{
    collections::HashMap,
    io::{BufRead, BufReader, Write},
    process::exit,
};

use regex::Regex;

mod symm;

use symm::Atom;

#[derive(Debug)]
struct Spectro {
    atoms: Vec<Atom>,
    freqs: Vec<f64>,
    disps: Vec<Vec<f64>>,
}

impl Spectro {
    fn new() -> Self {
        Self {
            atoms: Vec::new(),
            freqs: Vec::new(),
            disps: Vec::new(),
        }
    }
}

fn parse_infile(filename: &str) -> Spectro {
    let atomic_weights = HashMap::from([("1.0078250", 1), ("12.0000000", 6)]);
    let f = match std::fs::File::open(filename) {
        Ok(it) => it,
        Err(err) => {
            eprintln!("failed to open '{filename}' for parsing with '{err}'");
            exit(1);
        }
    };
    let lines = BufReader::new(f).lines().flatten();
    let mut skip = 0;
    let mut in_geom = false;
    let mut in_lxm = false;
    // the block of the LXM matrix
    let mut block = 0;
    let mut spectro = Spectro::new();
    let disp = Regex::new(r"^\d+$").unwrap();
    let header = Regex::new(r"^(\s*\d+)+\s*$").unwrap();
    for line in lines {
        if skip > 0 {
            skip -= 1
        } else if line.contains("MOLECULAR PRINCIPAL GEOMETRY") {
            skip = 2;
            in_geom = true;
        } else if in_geom {
            let fields: Vec<_> = line.split_whitespace().collect();
            if fields.len() == 0 {
                in_geom = false;
            } else {
                let atomic_number = match atomic_weights.get(fields[4]) {
                    Some(a) => *a,
                    None => panic!(
                        "atom with weight {} not found, tell Brent!",
                        fields[4]
                    ),
                };
                if let [x, y, z] = fields[1..=3]
                    .iter()
                    .map(|s| s.parse::<f64>().unwrap())
                    .collect::<Vec<_>>()[..]
                {
                    spectro.atoms.push(Atom::new(atomic_number, x, y, z));
                }
            }
        } else if line.contains("LXM MATRIX") {
            skip = 2;
            in_lxm = true;
        } else if in_lxm {
            let fields: Vec<_> = line.split_whitespace().collect();
            if fields.len() == 0 {
                skip = 1;
            } else if header.is_match(&line) {
                block += 1;
                continue;
            } else if line.contains("--------") {
                continue;
            } else if line.contains("LX MATRIX") {
                in_lxm = false;
            } else if disp.is_match(fields[0]) {
                for (i, d) in fields[1..].iter().enumerate() {
                    let idx = 10 * block + i;
                    if spectro.disps.len() <= idx {
                        spectro.disps.resize(idx + 1, vec![]);
                    }
                    spectro.disps[idx].push(f64::from_str(d).unwrap());
                }
            } else {
                spectro.freqs.extend(
                    fields
                        .iter()
                        .filter_map(|s| {
                            let f = s.parse::<f64>().unwrap();
                            if f > 1.0 {
                                Some(f)
                            } else {
                                None
                            }
                        })
                        .collect::<Vec<_>>(),
                );
            }
        }
    }
    spectro.freqs.reverse();
    spectro.disps.reverse();
    spectro.disps =
        spectro.disps[spectro.disps.len() - spectro.freqs.len()..].to_vec();
    spectro
}

fn write_outfile(filename: &str, spectro: &Spectro) {
    let mut f = match std::fs::File::create(filename) {
        Ok(it) => it,
        Err(err) => {
            eprintln!("failed to open '{filename}' for writing with '{err}'");
            exit(1);
        }
    };
    make_outfile(&mut f, spectro);
}

fn make_outfile<W: Write>(w: &mut W, spectro: &Spectro) {
    write!(
        w,
        " Entering Gaussian System, Link 0=/usr/local/g16/g16

                         Standard orientation:
 ---------------------------------------------------------------------
 Center     Atomic      Atomic             Coordinates (Angstroms)
 Number     Number       Type             X           Y           Z
 ---------------------------------------------------------------------
"
    )
    .unwrap();
    for (i, atom) in spectro.atoms.iter().enumerate() {
        writeln!(
            w,
            "{:7}{:11}{:12}{:16.6}{:12.6}{:12.6}",
            i + 1,
            atom.atomic_number,
            0,
            atom.x,
            atom.y,
            atom.z
        )
        .unwrap();
    }
    write!(
	w,
	" ---------------------------------------------------------------------
    49 basis functions,    92 primitive gaussians,    49 cartesian basis functions
    10 alpha electrons       10 beta electrons
 **********************************************************************

            Population analysis using the SCF Density.

 **********************************************************************
 Harmonic frequencies (cm**-1), IR intensities (KM/Mole), Raman scattering
 activities (A**4/AMU), depolarization ratios for plane and unpolarized
 incident light, reduced masses (AMU), force constants (mDyne/A),
 and normal coordinates:
").unwrap();

    for (i, chunk) in spectro.freqs.chunks(3).enumerate() {
        for j in 0..chunk.len() {
            write!(w, "{:>23}", 3 * i + j + 1).unwrap();
        }
        writeln!(w,).unwrap();
        for _ in chunk {
            write!(w, "{:>23}", "A1").unwrap();
        }
        write!(w, "\n Frequencies --{:>12.4}", chunk[0]).unwrap();
        for freq in &chunk[1..] {
            write!(w, "{:>23.4}", freq).unwrap();
        }
        writeln!(w,).unwrap();
        writeln!(w, " Red. masses --      1.0000                 1.0000                 1.0000").unwrap();
        writeln!(w, " Frc consts  --      1.0000                 1.0000                 1.0000").unwrap();
        writeln!(w, " IR Inten    --      1.0000                 1.0000                 1.0000").unwrap();
        writeln!(w, "  Atom  AN      X      Y      Z        X      Y      Z        X      Y      Z").unwrap();
        for (a, atom) in spectro.atoms.iter().enumerate() {
            write!(w, "{:6}{:4}  ", a + 1, atom.atomic_number).unwrap();
            for j in 0..chunk.len() {
                for f in &spectro.disps[3 * i + j][3 * a..3 * a + 3] {
                    write!(w, "{:7.2}", f).unwrap();
                }
                write!(w, "  ").unwrap();
            }
            writeln!(w,).unwrap();
        }
    }
}

fn main() {
    let args: Vec<_> = std::env::args().collect();
    let infile = match args.get(1) {
        Some(infile) => infile,
        None => {
            eprintln!("{}: not enough arguments", args[0]);
            exit(1);
        }
    };
    let atoms = parse_infile(infile);
    let outfile = match args.get(2) {
        Some(outfile) => outfile,
        None => "scope.out",
    };
    write_outfile(outfile, &atoms);
    make_outfile(&mut std::io::stdout(), &atoms);
}
