use criterion::{black_box, criterion_group, criterion_main, Criterion};

use noether::units::f64;

use noether::boundaries::{Pbc, NoBounds, BoundaryConditions};

fn bench_minimg(c: &mut Criterion) {
    let a_pos = black_box([3.5 * f64::NM, 4.7 * f64::NM, 1.2 * f64::NM]);
    let b_pos = black_box([2.3 * f64::NM, 1.9 * f64::NM, 0.4 * f64::NM]);
    let bounds = Pbc::rhombic_dodecahedral_xyhex(5.0 * f64::NM);

    c.bench_function("minimg dist2", move |b| b.iter(|| {
        bounds.dist2(a_pos, b_pos)
    }));

    let bounds = Pbc::rhombic_dodecahedral_xyhex(5.0 * f64::NM);

    c.bench_function("minimg dist", move |b| b.iter(|| {
        bounds.dist(a_pos, b_pos)
    }));

    let bounds = NoBounds;

    c.bench_function("straight dist2", move |b| b.iter(|| {
        bounds.dist2(a_pos, b_pos)
    }));

    let bounds = NoBounds;

    c.bench_function("straight dist", move |b| b.iter(|| {
        bounds.dist(a_pos, b_pos)
    }));
}

criterion_group!(benches, bench_minimg);
criterion_main!(benches);
