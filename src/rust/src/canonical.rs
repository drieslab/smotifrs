//! Canonical motif labeling.
//!
//! Mirrors `R/motifs_igraph.R`'s approximate canonical labeling: sort
//! vertices by `(degree-in-motif, color)`, then concatenate the color
//! tuple. We pad degree as 2 digits so byte-lexicographic sort agrees
//! with numeric sort even when more than ~10 different degrees show
//! up (paranoia for v0.2 size-5+).
//!
//! Orbit ids: vertices that share the same `(degree, color)` after
//! sorting share an orbit; orbit_id is the rank of the unique tuple
//! in lex order.

/// Sort 2 elements ascending by byte-lex order.
#[inline]
pub fn sort2<T: Ord + Clone>(a: &T, b: &T) -> [T; 2] {
    if a <= b {
        [a.clone(), b.clone()]
    } else {
        [b.clone(), a.clone()]
    }
}

/// Sort 3 elements ascending by byte-lex order.
#[inline]
pub fn sort3<T: Ord + Clone>(a: &T, b: &T, c: &T) -> [T; 3] {
    let mut out = [a.clone(), b.clone(), c.clone()];
    out.sort();
    out
}

/// Sort 4 elements ascending. Direct sort is fast enough for k=4.
#[inline]
pub fn sort4<T: Ord + Clone>(a: &T, b: &T, c: &T, d: &T) -> [T; 4] {
    let mut out = [a.clone(), b.clone(), c.clone(), d.clone()];
    out.sort();
    out
}

/// Sort 4 (deg, color) pairs by the lex `(deg, color)` order, return
/// the colors in sorted order.
#[inline]
pub fn sort4_by_deg_color<'a>(
    d: [u8; 4],
    c: [&'a str; 4],
) -> [&'a str; 4] {
    // Pair each (deg, color) with its index, sort, project colors.
    let mut idx: [usize; 4] = [0, 1, 2, 3];
    idx.sort_by(|&i, &j| match d[i].cmp(&d[j]) {
        std::cmp::Ordering::Equal => c[i].cmp(c[j]),
        other => other,
    });
    [c[idx[0]], c[idx[1]], c[idx[2]], c[idx[3]]]
}

/// Compute orbit ids for a vector of (deg, color) pairs. Vertices
/// sharing a (deg, color) tuple share an orbit. Returned orbit_id is
/// 1-based rank in sorted unique-tuple order.
pub fn orbit_ids(deg: &[u8], color: &[&str]) -> Vec<i32> {
    let n = deg.len();
    let mut keys: Vec<(u8, &str)> = deg.iter().zip(color.iter())
        .map(|(d, c)| (*d, *c))
        .collect();
    let mut uniq = keys.clone();
    uniq.sort();
    uniq.dedup();
    let mut out = Vec::with_capacity(n);
    for k in &keys {
        let r = uniq.binary_search(k).unwrap_or(0);
        out.push((r + 1) as i32);
    }
    // suppress "borrow of moved" by keeping `keys` alive
    let _ = &mut keys;
    out
}
