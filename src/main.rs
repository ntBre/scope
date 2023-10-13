use std::io::Write;

use summarize::{Summary, Recompute};

fn make_outfile<W: Write>(w: &mut W, spectro: &Summary) {
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
    for (i, atom) in spectro.geom.atoms.iter().enumerate() {
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

    let mut irreps = spectro.irreps.chunks(3);
    for (i, chunk) in spectro.harm.chunks(3).enumerate() {
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
        for (a, atom) in spectro.geom.atoms.iter().enumerate() {
            write!(w, "{:6}{:4}  ", a + 1, atom.atomic_number).unwrap();
            for j in 0..chunk.len() {
                for f in &spectro.lxm[3 * i + j][3 * a..3 * a + 3] {
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
            eprintln!("usage: scope FILENAME");
            return;
        }
    };
    let spectro = Summary::new(infile, Recompute::No);
    make_outfile(&mut std::io::stdout(), &spectro);
}
