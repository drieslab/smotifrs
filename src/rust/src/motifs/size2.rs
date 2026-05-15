//! Size-2 motif enumeration: every undirected edge is one instance.
//!
//! Mirrors `R/motifs_igraph.R::.find_motifs_size2()`. With colored = TRUE,
//! the motif id is `"size2_edge_<lo>-<hi>"` where `(lo, hi)` are the two
//! endpoint colors sorted byte-lexicographically. With colored = FALSE,
//! a single `"size2_edge"` motif id with `count = nrow(edges)`.
//!
//! Instance ids are `"e_<i>"` where `i` runs over instances in the same
//! deterministic order the R kernel uses (`setorder(motif_id, source,
//! target)` then `.I`).

use std::collections::HashMap;

use crate::canonical::orbit_ids;
use crate::graph::SpatialGraphRs;
use crate::motifs::{
    catalog_from_counts, EdgeMeta, IncidenceTable, MotifResult,
};

pub fn find_size2(
    g: &SpatialGraphRs,
    colored: bool,
    anchored_on: Option<&[u32]>,
) -> MotifResult {
    // Re-emit canonical edge list: one (lo, hi) per packed key. Sort so
    // the deterministic iteration order matches R.
    let mut edges: Vec<(u32, u32)> = g.edge_keys.iter()
        .map(|k| {
            let lo = (*k >> 32) as u32;
            let hi = (*k & 0xFFFF_FFFF) as u32;
            (lo, hi)
        })
        .collect();

    if let Some(anchors) = anchored_on {
        let anchor_set: ahash::AHashSet<u32> = anchors.iter().copied().collect();
        edges.retain(|(s, t)| anchor_set.contains(s) || anchor_set.contains(t));
    }

    // Tag each edge with its (motif_id, color_tuple) before sorting by
    // (motif_id, source, target) so instance_ids increment in the same
    // order R would produce them.
    struct Tagged {
        s: u32,
        t: u32,
        ca: String,
        cb: String,
        motif_id: String,
        color_tuple: Option<String>,
    }
    let mut tagged: Vec<Tagged> = edges.iter().map(|(s, t)| {
        let ca = g.color_of(*s).to_string();
        let cb = g.color_of(*t).to_string();
        let (motif_id, color_tuple) = if colored {
            let (lo, hi) = if ca <= cb { (&ca, &cb) } else { (&cb, &ca) };
            let ct = format!("{}-{}", lo, hi);
            (format!("size2_edge_{}", ct), Some(ct))
        } else {
            (String::from("size2_edge"), None)
        };
        Tagged { s: *s, t: *t, ca, cb, motif_id, color_tuple }
    }).collect();
    tagged.sort_by(|a, b| (&a.motif_id, a.s, a.t).cmp(&(&b.motif_id, b.s, b.t)));

    let mut counts: HashMap<String, i32> = HashMap::new();
    let mut info: HashMap<String, (i32, String, Option<String>)> = HashMap::new();
    let mut incidence = IncidenceTable::default();
    let mut edge_meta = EdgeMeta::default();

    for (i, e) in tagged.iter().enumerate() {
        let inst = format!("e_{}", i + 1);
        *counts.entry(e.motif_id.clone()).or_insert(0) += 1;
        info.entry(e.motif_id.clone()).or_insert_with(|| {
            (2, "edge".to_string(), e.color_tuple.clone())
        });

        let s_id = g.cell_ids[e.s as usize].clone();
        let t_id = g.cell_ids[e.t as usize].clone();
        edge_meta.instance_id.push(inst.clone());
        edge_meta.v1.push(s_id.clone());
        edge_meta.v2.push(t_id.clone());

        // orbit_id: both vertices have deg = 1 in motif. With coloring,
        // same color → orbit 1 for both; different colors → 1/2 by lex.
        let degs: [u8; 2] = [1, 1];
        let orbits: Vec<i32> = if colored {
            let cs = [e.ca.as_str(), e.cb.as_str()];
            orbit_ids(&degs, &cs)
        } else {
            vec![1, 1]
        };

        incidence.motif_id.push(e.motif_id.clone());
        incidence.instance_id.push(inst.clone());
        incidence.cell_id.push(s_id);
        incidence.orbit_id.push(orbits[0]);

        incidence.motif_id.push(e.motif_id.clone());
        incidence.instance_id.push(inst);
        incidence.cell_id.push(t_id);
        incidence.orbit_id.push(orbits[1]);
    }

    MotifResult {
        catalog: catalog_from_counts(&counts, &info),
        incidence,
        edge_meta: Some(edge_meta),
        triangle_meta: None,
        wedge_meta: None,
        size4_meta: HashMap::new(),
    }
}
