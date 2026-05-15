//! Size-4 motif enumeration via 3-subgraph extension + classification.
//!
//! Mirrors `R/motifs_igraph.R::.find_motifs_size4()`:
//!   1. Enumerate every connected 3-vertex induced subgraph (triangles
//!      ∪ wedges), each as a sorted `[u32; 3]`.
//!   2. Extend each by every neighbor of its three vertices not in the
//!      set; sort the resulting `[u32; 4]` to canonical lex order.
//!   3. Dedupe via a `u128` packed-key hash set.
//!   4. For each unique 4-set: count induced edges (6 hash lookups),
//!      compute degree sequence, classify into one of
//!      `{path, claw, cycle, paw, diamond, K4}`.
//!   5. Color tuple = sort 4 vertices by `(deg, color)` lex order, then
//!      concatenate the colors with `"-"`.

use std::collections::HashMap;

use ahash::AHashSet;

use crate::canonical::{orbit_ids, sort4_by_deg_color};
use crate::graph::SpatialGraphRs;
use crate::motifs::size3::{enumerate_triangles, enumerate_wedges};
use crate::motifs::{
    catalog_from_counts, IncidenceTable, MotifResult, Size4Meta,
};

#[inline]
fn pack4(w: [u32; 4]) -> u128 {
    ((w[0] as u128) << 96)
        | ((w[1] as u128) << 64)
        | ((w[2] as u128) << 32)
        | (w[3] as u128)
}

#[inline]
fn sort4_u32(mut x: [u32; 4]) -> [u32; 4] {
    x.sort_unstable();
    x
}

pub fn find_size4(
    g: &SpatialGraphRs,
    colored: bool,
    anchored_on: Option<&[u32]>,
) -> MotifResult {
    // 1. All connected 3-subgraphs as sorted v1 < v2 < v3.
    let mut three_sets: Vec<[u32; 3]> = Vec::new();
    for t in enumerate_triangles(g) {
        let mut s = t;
        s.sort_unstable();
        three_sets.push(s);
    }
    for w in enumerate_wedges(g) {
        let mut s = w;
        s.sort_unstable();
        three_sets.push(s);
    }
    // The same 3-set can appear once from triangles and never from
    // wedges (or vice versa), so dedupe.
    three_sets.sort_unstable();
    three_sets.dedup();

    // 2 + 3. Expand by neighbors, sort, dedupe.
    let mut seen: AHashSet<u128> = AHashSet::new();
    let mut quads: Vec<[u32; 4]> = Vec::new();
    for t in &three_sets {
        for &v in t {
            for &d in g.neighbors_of(v) {
                if d == t[0] || d == t[1] || d == t[2] {
                    continue;
                }
                let q = sort4_u32([t[0], t[1], t[2], d]);
                let key = pack4(q);
                if seen.insert(key) {
                    quads.push(q);
                }
            }
        }
    }
    drop(seen);

    let anchor_set: Option<AHashSet<u32>> = anchored_on.map(|a| a.iter().copied().collect());

    let mut counts: HashMap<String, i32> = HashMap::new();
    let mut info: HashMap<String, (i32, String, Option<String>)> = HashMap::new();
    let mut incidence = IncidenceTable::default();
    let mut size4_meta: HashMap<String, Size4Meta> = HashMap::new();

    // 4. Classify and tag each quad.
    struct QuadTag {
        w: [u32; 4],
        d: [u8; 4],
        class: &'static str,
        motif_id: String,
        color_tuple: Option<String>,
    }
    let mut tagged: Vec<QuadTag> = Vec::with_capacity(quads.len());
    for q in &quads {
        if let Some(set) = &anchor_set {
            if !q.iter().any(|v| set.contains(v)) {
                continue;
            }
        }
        let e = [
            g.has_edge(q[0], q[1]),
            g.has_edge(q[0], q[2]),
            g.has_edge(q[0], q[3]),
            g.has_edge(q[1], q[2]),
            g.has_edge(q[1], q[3]),
            g.has_edge(q[2], q[3]),
        ];
        let n_edges: u8 = e.iter().map(|b| *b as u8).sum();
        let d = [
            (e[0] as u8) + (e[1] as u8) + (e[2] as u8),
            (e[0] as u8) + (e[3] as u8) + (e[4] as u8),
            (e[1] as u8) + (e[3] as u8) + (e[5] as u8),
            (e[2] as u8) + (e[4] as u8) + (e[5] as u8),
        ];
        let mut ds = d;
        ds.sort_unstable();
        let class = match (n_edges, ds) {
            (3, [1, 1, 1, 3]) => "claw",
            (3, [1, 1, 2, 2]) => "path",
            (4, [2, 2, 2, 2]) => "cycle",
            (4, [1, 2, 2, 3]) => "paw",
            (5, [2, 2, 3, 3]) => "diamond",
            (6, [3, 3, 3, 3]) => "K4",
            _ => continue, // shouldn't happen for connected 4-subgraphs
        };

        let c0 = g.color_of(q[0]);
        let c1 = g.color_of(q[1]);
        let c2 = g.color_of(q[2]);
        let c3 = g.color_of(q[3]);
        let (motif_id, color_tuple) = if colored {
            let sorted = sort4_by_deg_color(d, [c0, c1, c2, c3]);
            let ct = format!("{}-{}-{}-{}", sorted[0], sorted[1], sorted[2], sorted[3]);
            (format!("size4_{}_{}", class, ct), Some(ct))
        } else {
            (format!("size4_{}", class), None)
        };

        tagged.push(QuadTag {
            w: *q, d, class, motif_id, color_tuple,
        });
    }
    tagged.sort_by(|a, b| {
        (&a.motif_id, a.w[0], a.w[1], a.w[2], a.w[3])
            .cmp(&(&b.motif_id, b.w[0], b.w[1], b.w[2], b.w[3]))
    });

    for (i, t) in tagged.iter().enumerate() {
        let inst = format!("q_{}", i + 1);
        *counts.entry(t.motif_id.clone()).or_insert(0) += 1;
        info.entry(t.motif_id.clone()).or_insert_with(|| {
            (4, t.class.to_string(), t.color_tuple.clone())
        });

        let cell = |idx: u32| g.cell_ids[idx as usize].clone();
        let ids = [cell(t.w[0]), cell(t.w[1]), cell(t.w[2]), cell(t.w[3])];

        // instance_meta per class.
        let bucket = size4_meta.entry(t.class.to_string()).or_default();
        bucket.instance_id.push(inst.clone());
        bucket.w1.push(ids[0].clone());
        bucket.w2.push(ids[1].clone());
        bucket.w3.push(ids[2].clone());
        bucket.w4.push(ids[3].clone());
        bucket.d1.push(t.d[0] as i32);
        bucket.d2.push(t.d[1] as i32);
        bucket.d3.push(t.d[2] as i32);
        bucket.d4.push(t.d[3] as i32);

        // orbit_id per vertex.
        let orbits: Vec<i32> = if colored {
            let colors = [
                g.color_of(t.w[0]),
                g.color_of(t.w[1]),
                g.color_of(t.w[2]),
                g.color_of(t.w[3]),
            ];
            orbit_ids(&t.d, &colors)
        } else {
            // Without color, orbit collapses to deg rank.
            let degs = [t.d[0], t.d[1], t.d[2], t.d[3]];
            let mut uniq = degs.to_vec();
            uniq.sort_unstable();
            uniq.dedup();
            degs.iter()
                .map(|d| (uniq.iter().position(|x| x == d).unwrap() + 1) as i32)
                .collect()
        };

        for k in 0..4 {
            incidence.motif_id.push(t.motif_id.clone());
            incidence.instance_id.push(inst.clone());
            incidence.cell_id.push(ids[k].clone());
            incidence.orbit_id.push(orbits[k]);
        }
    }

    MotifResult {
        catalog: catalog_from_counts(&counts, &info),
        incidence,
        edge_meta: None,
        triangle_meta: None,
        wedge_meta: None,
        size4_meta,
    }
}
