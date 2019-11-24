use crate::units::length::angstrom;
use crate::units::f64::*;

use chemfiles::{Trajectory, Frame};

use std::path::Path;
use std::fs::File;

use std::io;
use std::io::BufRead;

/// Creates a custom IO error with a message
fn custom_error<T>(message: &str) -> Result<T, io::Error> {
    Err(io::Error::new(io::ErrorKind::Other, message))
}

/// Reads a column from an xvg file on disk into a vector
///
/// Does not perform any numeric interpretation.
pub fn read_xvg<P: AsRef<Path>>(path: P, column: usize) -> Result<Vec<String>, io::Error> {
    let file = File::open(path)?;
    let reader = io::BufReader::new(file);
    let mut out = Vec::new();

    // Read the file line by line using the lines() iterator from std::io::BufRead.
    for line in reader.lines() {
        let line = line?;

        match line.trim().get(0..1) {
            None | Some("@") | Some("#") => continue,
            Some(_) => ()
        };

        match line.split_whitespace().nth(column) {
            None => return custom_error("A row has too few values"),
            Some(data) => out.push(data.to_string())
        }
    }

    Ok(out)
}

/// Reads positions into a vector from a file in any format supported by chemfile
pub fn read_positions<P: AsRef<Path>>(path: P) -> Result<Vec<Positions>, io::Error> {
    let mut pos = Vec::new();
    let mut frame = Frame::new();

    let mut traj = match Trajectory::open(path, 'r') {
        Err(_) => return custom_error("File could not be opened"),
        Ok(traj) => traj,
    };

    loop {
        match traj.read(&mut frame) {
            Ok(()) => (),
            Err(_) => break,
        }

        pos.push(
            frame
                .positions()
                .iter()
                .map(|x| [
                    Length::new::<angstrom>(x[0]),
                    Length::new::<angstrom>(x[1]),
                    Length::new::<angstrom>(x[2])
                ])
                .collect()
        )
    }

    Ok(pos)
}
