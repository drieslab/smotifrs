//! Exact canonical forms for colored motifs of size 2, 3 and 4.
//!
//! The previous scheme sorted vertices by `(degree, color)`, which is not a
//! canonical form: a colored C4 `A-B-A-B` and `A-A-B-B` both reduce to degrees
//! `(2,2,2,2)` and sorted colors `A-A-B-B`, so two structurally different
//! motifs collided on one id. The same happens for diamonds and paws.
//!
//! Here the canonical form is exact, and nothing is hard-coded per topology:
//!
//! * A k-subset's structure is a bitmask over its `k*(k-1)/2` vertex pairs.
//! * The **minimum** of that mask over all `k!` vertex permutations is a
//!   canonical form for the unlabeled graph, and the argmin gives a slot
//!   ordering. Two subsets are isomorphic iff their canonical masks are equal.
//! * The automorphism group of a topology is exactly the set of permutations
//!   that fix its canonical mask. It is derived, never tabulated by hand.
//! * A colored motif's canonical color key is the lexicographic minimum of the
//!   slot-ordered color tuple over that automorphism group.
//!
//! Costs are split deliberately. Structure canonicalization is the expensive
//! part (up to 24 permutations) and runs **once per instance**, at enumeration
//! time. The color key runs once per instance *per permutation draw*, and only
//! minimizes over `Aut`, which has at most 24 and typically 2-8 elements.
//!
//! Colors are integer codes throughout; no strings are built in any hot loop.
//! Packing four 16-bit codes big-endian into a `u64` makes numeric order on the
//! packed value agree with lexicographic order on the tuple.

/// Permutations of `0..k` for k = 2, 3, 4. Indexed by `PERMS[k]`.
const PERMS2: [[u8; 4]; 2] = [[0, 1, 2, 3], [1, 0, 2, 3]];

const PERMS3: [[u8; 4]; 6] = [
    [0, 1, 2, 3], [0, 2, 1, 3], [1, 0, 2, 3],
    [1, 2, 0, 3], [2, 0, 1, 3], [2, 1, 0, 3],
];

const PERMS4: [[u8; 4]; 24] = [
    [0, 1, 2, 3], [0, 1, 3, 2], [0, 2, 1, 3], [0, 2, 3, 1],
    [0, 3, 1, 2], [0, 3, 2, 1], [1, 0, 2, 3], [1, 0, 3, 2],
    [1, 2, 0, 3], [1, 2, 3, 0], [1, 3, 0, 2], [1, 3, 2, 0],
    [2, 0, 1, 3], [2, 0, 3, 1], [2, 1, 0, 3], [2, 1, 3, 0],
    [2, 3, 0, 1], [2, 3, 1, 0], [3, 0, 1, 2], [3, 0, 2, 1],
    [3, 1, 0, 2], [3, 1, 2, 0], [3, 2, 0, 1], [3, 2, 1, 0],
];

#[inline]
pub fn perms_for(k: usize) -> &'static [[u8; 4]] {
    match k {
        2 => &PERMS2,
        3 => &PERMS3,
        4 => &PERMS4,
        _ => panic!("unsupported motif size {}", k),
    }
}

/// Bit position of the unordered pair `(i, j)`, `i < j`, for motif size `k`.
/// Layout for k=4: (0,1)=0 (0,2)=1 (0,3)=2 (1,2)=3 (1,3)=4 (2,3)=5.
#[inline]
pub fn pair_bit(i: usize, j: usize, k: usize) -> usize {
    let (a, b) = if i < j { (i, j) } else { (j, i) };
    let mut off = 0usize;
    for r in 0..a {
        off += k - 1 - r;
    }
    off + (b - a - 1)
}

/// Rewrite an adjacency mask under a vertex permutation: the returned mask has
/// the bit for `(i, j)` set iff the input had an edge between `perm[i]` and
/// `perm[j]`.
#[inline]
pub fn perm_adj(adj: u8, perm: &[u8; 4], k: usize) -> u8 {
    let mut out = 0u8;
    for i in 0..k {
        for j in (i + 1)..k {
            let src = pair_bit(perm[i] as usize, perm[j] as usize, k);
            if adj & (1u8 << src) != 0 {
                out |= 1u8 << pair_bit(i, j, k);
            }
        }
    }
    out
}

/// Canonical structure of a k-subset: the lexicographically smallest adjacency
/// mask over all vertex permutations, together with the permutation achieving
/// it (`slots[i]` = which input vertex position sits in canonical slot `i`).
#[inline]
pub fn canonical_structure(adj: u8, k: usize) -> (u8, [u8; 4]) {
    let mut best = u8::MAX;
    let mut best_perm = [0u8, 1, 2, 3];
    for p in perms_for(k) {
        let m = perm_adj(adj, p, k);
        if m < best {
            best = m;
            best_perm = *p;
        }
    }
    (best, best_perm)
}

/// Automorphism group of the topology whose canonical mask is `canon`: every
/// permutation that leaves the mask unchanged. Derived, not tabulated.
pub fn automorphisms(canon: u8, k: usize) -> Vec<[u8; 4]> {
    perms_for(k)
        .iter()
        .filter(|p| perm_adj(canon, p, k) == canon)
        .copied()
        .collect()
}

/// Pack up to four 16-bit color codes big-endian, so numeric order on the
/// result matches lexicographic order on the tuple.
#[inline]
pub fn pack_colors(c: &[u32; 4], k: usize) -> u64 {
    let mut out = 0u64;
    for i in 0..k {
        out |= (c[i] as u64 & 0xFFFF) << (48 - 16 * i);
    }
    out
}

#[inline]
pub fn unpack_colors(key: u64, k: usize) -> [u32; 4] {
    let mut out = [0u32; 4];
    for i in 0..k {
        out[i] = ((key >> (48 - 16 * i)) & 0xFFFF) as u32;
    }
    out
}

/// Canonical color key for one instance: the smallest packed color tuple over
/// the topology's automorphism group. `slot_colors` must already be in
/// canonical slot order.
#[inline]
pub fn canonical_color_key(slot_colors: &[u32; 4], auts: &[[u8; 4]], k: usize) -> u64 {
    let mut best = u64::MAX;
    for a in auts {
        let mut c = [0u32; 4];
        for i in 0..k {
            c[i] = slot_colors[a[i] as usize];
        }
        let p = pack_colors(&c, k);
        if p < best {
            best = p;
        }
    }
    best
}

/// The nine connected topologies on 2, 3 and 4 vertices, discovered by
/// enumerating every adjacency mask and keeping the connected canonical ones.
#[derive(Debug, Clone)]
pub struct Topology {
    pub canon: u8,
    pub k: usize,
    pub auts: Vec<[u8; 4]>,
    pub name: &'static str,
}

fn is_connected(adj: u8, k: usize) -> bool {
    let mut seen = [false; 4];
    let mut stack = vec![0usize];
    seen[0] = true;
    let mut n = 1;
    while let Some(v) = stack.pop() {
        for u in 0..k {
            if u != v && !seen[u] && adj & (1u8 << pair_bit(v, u, k)) != 0 {
                seen[u] = true;
                n += 1;
                stack.push(u);
            }
        }
    }
    n == k
}

/// Name a size-4 topology from its edge count and sorted degree sequence, so
/// the emitted ids stay compatible with the existing `size4_<class>` grammar.
fn name_of(canon: u8, k: usize) -> &'static str {
    let mut deg = [0u8; 4];
    for i in 0..k {
        for j in 0..k {
            if i != j && canon & (1u8 << pair_bit(i, j, k)) != 0 {
                deg[i] += 1;
            }
        }
    }
    let mut d = deg[..k].to_vec();
    d.sort_unstable();
    let ne = (canon.count_ones()) as usize;
    match (k, ne, d.as_slice()) {
        (2, 1, _) => "edge",
        (3, 2, _) => "open",
        (3, 3, _) => "closed",
        (4, 3, [1, 1, 1, 3]) => "claw",
        (4, 3, [1, 1, 2, 2]) => "path",
        (4, 4, [2, 2, 2, 2]) => "cycle",
        (4, 4, [1, 2, 2, 3]) => "paw",
        (4, 5, _) => "diamond",
        (4, 6, _) => "K4",
        _ => "unknown",
    }
}

/// All connected topologies for a given motif size, keyed by canonical mask.
pub fn topologies(k: usize) -> Vec<Topology> {
    let nbits = k * (k - 1) / 2;
    let mut seen: Vec<u8> = Vec::new();
    let mut out = Vec::new();
    for adj in 0u16..(1u16 << nbits) {
        let adj = adj as u8;
        if !is_connected(adj, k) {
            continue;
        }
        let (canon, _) = canonical_structure(adj, k);
        if seen.contains(&canon) {
            continue;
        }
        seen.push(canon);
        out.push(Topology {
            canon,
            k,
            auts: automorphisms(canon, k),
            name: name_of(canon, k),
        });
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pair_bits_are_a_bijection() {
        for k in 2..=4 {
            let mut bits = Vec::new();
            for i in 0..k {
                for j in (i + 1)..k {
                    bits.push(pair_bit(i, j, k));
                }
            }
            let mut s = bits.clone();
            s.sort_unstable();
            s.dedup();
            assert_eq!(s.len(), k * (k - 1) / 2, "k={}", k);
            assert_eq!(*s.iter().max().unwrap(), k * (k - 1) / 2 - 1);
        }
    }

    #[test]
    fn there_are_nine_connected_topologies() {
        assert_eq!(topologies(2).len(), 1);
        assert_eq!(topologies(3).len(), 2);
        assert_eq!(topologies(4).len(), 6);
    }

    #[test]
    fn automorphism_group_orders_are_textbook() {
        let by = |k: usize, n: &str| -> usize {
            topologies(k).into_iter().find(|t| t.name == n).unwrap().auts.len()
        };
        assert_eq!(by(2, "edge"), 2); // S2
        assert_eq!(by(3, "open"), 2); // reverse the wedge
        assert_eq!(by(3, "closed"), 6); // S3
        assert_eq!(by(4, "claw"), 6); // S3 on the leaves
        assert_eq!(by(4, "path"), 2); // reverse
        assert_eq!(by(4, "cycle"), 8); // D4
        assert_eq!(by(4, "paw"), 2);
        assert_eq!(by(4, "diamond"), 4);
        assert_eq!(by(4, "K4"), 24); // S4
    }

    #[test]
    fn isomorphic_subsets_share_a_canonical_mask() {
        // every relabeling of a paw must canonicalize identically
        let paw: u8 = (1 << pair_bit(0, 1, 4))
            | (1 << pair_bit(1, 2, 4))
            | (1 << pair_bit(0, 2, 4))
            | (1 << pair_bit(0, 3, 4));
        let (c0, _) = canonical_structure(paw, 4);
        for p in perms_for(4) {
            let (c, _) = canonical_structure(perm_adj(paw, p, 4), 4);
            assert_eq!(c, c0);
        }
    }

    // The regression this module exists for: on a 4-cycle, the colorings
    // A-B-A-B (alternating) and A-A-B-B (adjacent pairs) are NOT isomorphic as
    // colored graphs. The old (degree, color) sort gave both the same id.
    #[test]
    fn alternating_and_adjacent_c4_colorings_differ() {
        let c4: u8 = (1 << pair_bit(0, 1, 4))
            | (1 << pair_bit(1, 2, 4))
            | (1 << pair_bit(2, 3, 4))
            | (1 << pair_bit(0, 3, 4));
        let (canon, slots) = canonical_structure(c4, 4);
        let auts = automorphisms(canon, 4);

        let key = |cols: [u32; 4]| {
            let mut sc = [0u32; 4];
            for i in 0..4 {
                sc[i] = cols[slots[i] as usize];
            }
            canonical_color_key(&sc, &auts, 4)
        };

        let alternating = key([0, 1, 0, 1]); // A B A B around the cycle
        let adjacent = key([0, 0, 1, 1]); // A A B B around the cycle
        assert_ne!(
            alternating, adjacent,
            "alternating and adjacent C4 colorings must not collide"
        );

        // ...while genuine rotations/reflections of one coloring must agree
        assert_eq!(alternating, key([1, 0, 1, 0]));
        assert_eq!(adjacent, key([0, 1, 1, 0]));
        assert_eq!(adjacent, key([1, 1, 0, 0]));
    }

    #[test]
    fn diamond_colorings_respect_degree_roles() {
        let dia: u8 = (1 << pair_bit(0, 1, 4))
            | (1 << pair_bit(0, 2, 4))
            | (1 << pair_bit(0, 3, 4))
            | (1 << pair_bit(1, 2, 4))
            | (1 << pair_bit(1, 3, 4));
        let (canon, slots) = canonical_structure(dia, 4);
        let auts = automorphisms(canon, 4);
        let key = |cols: [u32; 4]| {
            let mut sc = [0u32; 4];
            for i in 0..4 {
                sc[i] = cols[slots[i] as usize];
            }
            canonical_color_key(&sc, &auts, 4)
        };
        // colouring the two degree-3 vertices differs from colouring the two
        // degree-2 vertices, even though the color multiset is identical
        assert_ne!(key([0, 0, 1, 1]), key([1, 1, 0, 0]));
        // swapping within a degree class is an automorphism
        assert_eq!(key([0, 1, 2, 3]), key([1, 0, 2, 3]));
        assert_eq!(key([0, 1, 2, 3]), key([0, 1, 3, 2]));
    }

    #[test]
    fn packing_preserves_lexicographic_order() {
        assert!(pack_colors(&[0, 0, 0, 1], 4) < pack_colors(&[0, 0, 1, 0], 4));
        assert!(pack_colors(&[1, 0, 0, 0], 4) > pack_colors(&[0, 9, 9, 9], 4));
        let c = [3u32, 1, 4, 1];
        assert_eq!(unpack_colors(pack_colors(&c, 4), 4), c);
    }
}
