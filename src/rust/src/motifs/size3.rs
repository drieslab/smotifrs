//! Size-3 motif enumeration: closed (triangles) and open (wedges).
//!
//! - Triangles: sorted-merge intersection of neighbor lists. Each
//!   triple (u, v, w) is emitted once with `u < v < w`.
//! - Wedges: for each center b, every unordered pair (a, c) of
//!   neighbors with a < c such that (a, c) is **not** an edge.
//!
//! Output schema must match `R/motifs_igraph.R`'s
//! `.find_motifs_size3_closed()` and `.find_motifs_size3_open()` —
//! same motif_id strings, same instance_id prefixes (`tri_`, `wedge_`),
//! same role-canonicalized instance_meta (triangle: lex-sorted; wedge:
//! end1 < end2 by cell-id index, center separate).

use std::collections::HashMap;

use crate::canonical::orbit_ids;
use crate::graph::SpatialGraphRs;
use crate::motifs::{
    catalog_from_counts, IncidenceTable, MotifResult, TriangleMeta, WedgeMeta,
};

pub fn find_size3(
    g: &SpatialGraphRs,
    colored: bool,
    anchored_on: Option<&[u32]>,
) -> MotifResult {
    let triangles = enumerate_triangles(g);
    let wedges = enumerate_wedges(g);

    let anchor_set: Option<ahash::AHashSet<u32>> = anchored_on.map(|a| {
        a.iter().copied().collect()
    });

    let mut counts: HashMap<String, i32> = HashMap::new();
    let mut info: HashMap<String, (i32, String, Option<String>)> = HashMap::new();
    let mut incidence = IncidenceTable::default();
    let mut tri_meta = TriangleMeta::default();
    let mut wedge_meta = WedgeMeta::default();

    // ---- triangles ------------------------------------------------------

    // Tag and sort by (motif_id, v1, v2, v3) for deterministic instance ids.
    struct TriTag {
        v: [u32; 3],
        c: [String; 3],
        motif_id: String,
        color_tuple: Option<String>,
    }
    let mut tri_tagged: Vec<TriTag> = triangles.iter().filter_map(|tri| {
        if let Some(set) = &anchor_set {
            if !tri.iter().any(|v| set.contains(v)) {
                return None;
            }
        }
        let c1 = g.color_of(tri[0]).to_string();
        let c2 = g.color_of(tri[1]).to_string();
        let c3 = g.color_of(tri[2]).to_string();
        let (motif_id, color_tuple) = if colored {
            let mut sorted = [c1.clone(), c2.clone(), c3.clone()];
            sorted.sort();
            let ct = format!("{}-{}-{}", sorted[0], sorted[1], sorted[2]);
            (format!("size3_closed_{}", ct), Some(ct))
        } else {
            (String::from("size3_closed"), None)
        };
        Some(TriTag { v: *tri, c: [c1, c2, c3], motif_id, color_tuple })
    }).collect();
    tri_tagged.sort_by(|a, b| {
        (&a.motif_id, a.v[0], a.v[1], a.v[2])
            .cmp(&(&b.motif_id, b.v[0], b.v[1], b.v[2]))
    });

    for (i, t) in tri_tagged.iter().enumerate() {
        let inst = format!("tri_{}", i + 1);
        *counts.entry(t.motif_id.clone()).or_insert(0) += 1;
        info.entry(t.motif_id.clone()).or_insert_with(|| {
            (3, "triangle".to_string(), t.color_tuple.clone())
        });

        let v1 = g.cell_ids[t.v[0] as usize].clone();
        let v2 = g.cell_ids[t.v[1] as usize].clone();
        let v3 = g.cell_ids[t.v[2] as usize].clone();
        tri_meta.instance_id.push(inst.clone());
        tri_meta.v1.push(v1.clone());
        tri_meta.v2.push(v2.clone());
        tri_meta.v3.push(v3.clone());

        // orbit_id: all vertices have deg = 2.
        let degs: [u8; 3] = [2, 2, 2];
        let orbits: Vec<i32> = if colored {
            let cs = [t.c[0].as_str(), t.c[1].as_str(), t.c[2].as_str()];
            orbit_ids(&degs, &cs)
        } else {
            vec![1, 1, 1]
        };

        for k in 0..3 {
            incidence.motif_id.push(t.motif_id.clone());
            incidence.instance_id.push(inst.clone());
            incidence.cell_id.push([&v1, &v2, &v3][k].clone());
            incidence.orbit_id.push(orbits[k]);
        }
    }

    // ---- wedges ---------------------------------------------------------

    struct WedgeTag {
        end_a: u32,
        center: u32,
        end_c: u32,
        ca: String,
        cb: String,
        cc: String,
        motif_id: String,
        color_tuple: Option<String>,
    }
    let mut wedge_tagged: Vec<WedgeTag> = wedges.iter().filter_map(|w| {
        if let Some(set) = &anchor_set {
            if !(set.contains(&w[0]) || set.contains(&w[1]) || set.contains(&w[2])) {
                return None;
            }
        }
        let ca = g.color_of(w[0]).to_string();
        let cb = g.color_of(w[1]).to_string();
        let cc = g.color_of(w[2]).to_string();
        let (motif_id, color_tuple) = if colored {
            let (lo, hi) = if ca <= cc { (ca.clone(), cc.clone()) }
                           else { (cc.clone(), ca.clone()) };
            let ct = format!("{}-{}_{}", lo, hi, cb);
            (format!("size3_open_{}", ct), Some(ct))
        } else {
            (String::from("size3_open"), None)
        };
        Some(WedgeTag {
            end_a: w[0], center: w[1], end_c: w[2],
            ca, cb, cc, motif_id, color_tuple,
        })
    }).collect();
    wedge_tagged.sort_by(|a, b| {
        (&a.motif_id, a.center, a.end_a, a.end_c)
            .cmp(&(&b.motif_id, b.center, b.end_a, b.end_c))
    });

    for (i, w) in wedge_tagged.iter().enumerate() {
        let inst = format!("wedge_{}", i + 1);
        *counts.entry(w.motif_id.clone()).or_insert(0) += 1;
        info.entry(w.motif_id.clone()).or_insert_with(|| {
            (3, "wedge".to_string(), w.color_tuple.clone())
        });

        let id_a = g.cell_ids[w.end_a as usize].clone();
        let id_b = g.cell_ids[w.center as usize].clone();
        let id_c = g.cell_ids[w.end_c as usize].clone();
        wedge_meta.instance_id.push(inst.clone());
        wedge_meta.end1.push(id_a.clone());
        wedge_meta.end2.push(id_c.clone());
        wedge_meta.center.push(id_b.clone());

        // orbit_id: ends have deg = 1, center has deg = 2.
        let degs: [u8; 3] = [1, 2, 1];
        let orbits: Vec<i32> = if colored {
            let cs = [w.ca.as_str(), w.cb.as_str(), w.cc.as_str()];
            orbit_ids(&degs, &cs)
        } else {
            // R uses ranks of unique deg values: deg=1 → 1, deg=2 → 2.
            vec![1, 2, 1]
        };

        // Incidence emits in the canonical order R uses: end-a, center, end-c.
        let bundle = [(id_a, orbits[0]), (id_b, orbits[1]), (id_c, orbits[2])];
        for (cid, o) in bundle {
            incidence.motif_id.push(w.motif_id.clone());
            incidence.instance_id.push(inst.clone());
            incidence.cell_id.push(cid);
            incidence.orbit_id.push(o);
        }
    }

    MotifResult {
        catalog: catalog_from_counts(&counts, &info),
        incidence,
        edge_meta: None,
        triangle_meta: Some(tri_meta),
        wedge_meta: Some(wedge_meta),
        size4_meta: HashMap::new(),
    }
}

/// Sorted-merge triangle enumeration. Returns each triangle as a
/// sorted `[u32; 3]` with `v1 < v2 < v3`. Linear-scan intersections
/// take advantage of pre-sorted neighbor slices in the CSR.
pub fn enumerate_triangles(g: &SpatialGraphRs) -> Vec<[u32; 3]> {
    let mut out: Vec<[u32; 3]> = Vec::new();
    for u in 0..g.n_nodes as u32 {
        let n_u = g.neighbors_of(u);
        // Iterate v ∈ N(u) with v > u so each triangle is found once at
        // its smallest vertex.
        for &v in n_u.iter().filter(|&&x| x > u) {
            let n_v = g.neighbors_of(v);
            // Two-pointer intersection of N(u) and N(v), keeping
            // elements > v so triangle (u, v, w) has u < v < w.
            let mut i = n_u.iter().position(|&x| x > v).unwrap_or(n_u.len());
            let mut j = n_v.iter().position(|&x| x > v).unwrap_or(n_v.len());
            while i < n_u.len() && j < n_v.len() {
                let a = n_u[i];
                let b = n_v[j];
                if a == b {
                    out.push([u, v, a]);
                    i += 1;
                    j += 1;
                } else if a < b {
                    i += 1;
                } else {
                    j += 1;
                }
            }
        }
    }
    out
}

/// Wedge enumeration: every (a, b, c) where b is a common neighbor of
/// a and c, a < c (canonicalize), and (a, c) is **not** an edge (else
/// it'd be a triangle).
pub fn enumerate_wedges(g: &SpatialGraphRs) -> Vec<[u32; 3]> {
    let mut out: Vec<[u32; 3]> = Vec::new();
    for b in 0..g.n_nodes as u32 {
        let n_b = g.neighbors_of(b);
        for i in 0..n_b.len() {
            let a = n_b[i];
            for j in (i + 1)..n_b.len() {
                let c = n_b[j];
                // n_b is sorted, so a < c by construction. We need
                // a != b and c != b — they are, since neighbor lists
                // exclude self in undirected canonical edges.
                if !g.has_edge(a, c) {
                    out.push([a, b, c]);
                }
            }
        }
    }
    out
}
