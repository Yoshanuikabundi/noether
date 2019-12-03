use crate::units::velocity::angstrom_per_picosecond;
use crate::units::length::angstrom;
use crate::units::f64::*;
use crate::units::Quantity;

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
/// I assume here that chemfiles reliably produces angstrom units
pub fn read_positions<P: AsRef<Path>>(path: P) -> Result<Vec<Positions>, io::Error> {
    read_frames(path, &Frame::positions, &Length::new::<angstrom>)
}

/// Reads velocities into a vector from a file in any format supported by chemfile
/// I assume here that chemfiles reliably produces angstrom per picosecond units
pub fn read_velocities<P: AsRef<Path>>(path: P) -> Result<Vec<Velocities>, io::Error> {
    read_frames(path, &Frame::velocities, &Velocity::new::<angstrom_per_picosecond>)
}

/// Read a property from a frame into a vector with units
fn read_frame<D: ?Sized, U: ?Sized, V>(
    frame: &Frame,
    prop_func: &dyn Fn(&chemfiles::Frame) -> &[[V; 3]],
    unit_constructor: &dyn Fn(V) -> Quantity<D, U, V>
    ) -> Vec<[Quantity<D, U, V>; 3]>
    where
        D: crate::units::Dimension,
        U: crate::units::Units<V>,
        V: uom::num_traits::Num + uom::Conversion<V> + Copy
{
    prop_func(&frame)
        .iter()
        .map(|&[x, y, z]| [
            unit_constructor(x),
            unit_constructor(y),
            unit_constructor(z)
        ])
        .collect()
}

/// Reads some property of all frames into a vector from a trajectory file in any format supported by chemfile
fn read_frames<P: AsRef<Path>, D: ?Sized, U: ?Sized, V>(
    path: P,
    prop_func: &dyn Fn(&Frame) -> &[[V; 3]],
    unit_constructor: &dyn Fn(V) -> Quantity<D, U, V>
    ) -> Result<Vec<Vec<[Quantity<D, U, V>; 3]>>, io::Error>
    where
        D: crate::units::Dimension,
        U: crate::units::Units<V>,
        V: uom::num_traits::Num + uom::Conversion<V> + Copy
{
    // Open the trajectory file
    let mut traj = match Trajectory::open(path, 'r') {
        Err(_) => return custom_error("File could not be opened"),
        Ok(traj) => traj,
    };

    // Create a frame object to read with
    let mut frame = Frame::new();
    match traj.nsteps() {
        // If trajectory has a known number of steps, iterate over them
        Ok(n) => {
            let mut pos = Vec::with_capacity(n as usize);
            for _ in 0..n {
                match traj.read(&mut frame) { // Only this function has access to traj, so we can safely use read
                    Ok(()) => (),
                    Err(_) => return custom_error("Read failed in middle of trajectory"),
                };

                pos.push(read_frame(&frame, prop_func, unit_constructor));
            }
            Ok(pos)
        },
        // Otherwise, iteratively read from the trajectory until it errs
        Err(_) => {
            let mut pos = Vec::new();

            loop {
                match traj.read(&mut frame) {
                    Ok(()) => (),
                    Err(_) => break,
                };

                pos.push(read_frame(&frame, prop_func, unit_constructor));
            }
            Ok(pos)
        }
    }

}
