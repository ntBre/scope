use std::{
    collections::HashMap,
    io::{BufRead, BufReader, Write},
    process::exit,
};

#[derive(Debug)]
struct Atom(usize, f64, f64, f64);

fn parse_infile(filename: &str) -> Vec<Atom> {
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
    let mut atoms = Vec::new();
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
                    atoms.push(Atom(atomic_number, x, y, z));
                }
            }
        }
    }
    atoms
}

fn write_outfile(filename: &str, atoms: &Vec<Atom>) {
    let mut f = match std::fs::File::create(filename) {
        Ok(it) => it,
        Err(err) => {
            eprintln!("failed to open '{filename}' for writing with '{err}'");
            exit(1);
        }
    };
    make_outfile(&mut f, atoms);
}

fn make_outfile<W: Write>(w: &mut W, atoms: &Vec<Atom>) {
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
    for (i, atom) in atoms.iter().enumerate() {
        writeln!(
            w,
            "{:7}{:11}{:12}{:16.6}{:12.6}{:12.6}",
            i + 1,
            atom.0,
            0,
            atom.1,
            atom.2,
            atom.3
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
    // make_outfile(&mut std::io::stdout(), &atoms);
}
