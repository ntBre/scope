use std::str::FromStr;
use std::{
    collections::HashMap,
    io::{BufRead, BufReader, Write},
    process::exit,
};

use regex::Regex;

mod symm;

use symm::{Atom, Irrep, Molecule};

/// threshold for discarding rotations and translations
const ROTRANS_THRSH: f64 = 30.0;

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

    fn load(filename: &str) -> Self {
        let atomic_weights = HashMap::from([
            ("1.0078250", 1),
            ("4.0026032", 2),
            ("7.0160030", 3),
            ("9.0121822", 4),
            ("11.0093054", 5),
            ("12.0000000", 6),
            ("14.0030740", 7),
            ("15.9949146", 8),
            ("18.9984032", 9),
            ("19.9924356", 10),
            ("22.9897677", 11),
            ("23.9850423", 12),
            ("26.9815386", 13),
            ("27.9769271", 14),
            ("30.9737620", 15),
            ("31.9720707", 16),
            ("34.9688527", 17),
            ("39.9623837", 18),
        ]);
        let f = match std::fs::File::open(filename) {
            Ok(it) => it,
            Err(err) => {
                eprintln!(
                    "failed to open '{filename}' for parsing with '{err}'"
                );
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
                // reset these. for degmodes it gets printed twice
                block = 0;
                spectro.disps = Vec::new();
                spectro.freqs = Vec::new();
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
                                if f > ROTRANS_THRSH {
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
}

fn make_outfile<W: Write>(w: &mut W, spectro: &Spectro, irreps: &Vec<Irrep>) {
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

    let mut irreps = irreps.chunks(3);
    for (i, chunk) in spectro.freqs.chunks(3).enumerate() {
        let irrep = irreps.next().unwrap();
        for j in 0..chunk.len() {
            write!(w, "{:>23}", 3 * i + j + 1).unwrap();
        }
        writeln!(w,).unwrap();
        for r in irrep {
            write!(w, "{:>23}", r.to_string()).unwrap();
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
    let spectro = Spectro::load(infile);
    let mol = Molecule::new(spectro.atoms.clone());
    let pg = mol.point_group();
    let mut irreps = Vec::new();
    for disp in &spectro.disps {
        let mol = mol.clone() + disp.clone();
        irreps.push(mol.irrep(&pg));
    }
    make_outfile(&mut std::io::stdout(), &spectro, &irreps);
}
