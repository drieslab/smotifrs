//! Compressed sparse-row spatial graph.
//!
//! `SpatialGraphRs` is the in-memory representation our motif kernels
//! work against. The build step canonicalizes the input edge list
//! (source < target, deduped, no self-loops) and produces, for every
//! vertex, a sorted slice of u32 neighbor indices. Sorted neighbor
//! lists let triangle/wedge/4-motif enumeration use linear-scan
//! intersections instead of hash lookups, which is what makes the
//! Rust backend competitive at million-cell scale.

use ahash::AHashSet;

#[derive(Debug)]
pub struct SpatialGraphRs {
    pub n_nodes: usize,
    pub cell_ids: Vec<String>,
    /// Neighbors of vertex `v` are `neighbors[neigh_offsets[v] .. neigh_offsets[v + 1]]`.
    pub neigh_offsets: Vec<u32>,
    pub neighbors: Vec<u32>,
    /// Per-vertex integer color code. Decode via `color_levels[c as usize]`.
    pub colors: Vec<u32>,
    pub color_levels: Vec<String>,
    /// Set of `(min_v, max_v)` packed into u64 for O(1) edge-membership lookup.
    /// Used by wedge/non-edge tests in size-3 open and size-4 enumeration.
    pub edge_keys: AHashSet<u64>,
}

#[inline]
pub fn pack_edge(a: u32, b: u32) -> u64 {
    let (lo, hi) = if a < b { (a, b) } else { (b, a) };
    ((lo as u64) << 32) | (hi as u64)
}

impl SpatialGraphRs {
    /// Build from raw vectors. `source`/`target` are u32 indices into
    /// `cell_ids` (already index-encoded). Edges are validated as
    /// `source < target` and deduped.
    pub fn build(
        cell_ids: Vec<String>,
        cell_types: Vec<String>,
        source_idx: &[u32],
        target_idx: &[u32],
    ) -> Result<Self, String> {
        let n = cell_ids.len();
        if cell_types.len() != n {
            return Err(format!(
                "cell_types length {} != cell_ids length {}",
                cell_types.len(),
                n
            ));
        }
        if source_idx.len() != target_idx.len() {
            return Err("source/target length mismatch".to_string());
        }

        // Encode cell_type strings into integer codes (sorted-unique levels).
        let mut levels: Vec<String> = cell_types.iter().cloned().collect();
        levels.sort();
        levels.dedup();
        let mut colors = Vec::with_capacity(n);
        for ct in &cell_types {
            let idx = levels.binary_search(ct).map_err(|_| {
                format!("internal: cell_type '{}' not in levels", ct)
            })?;
            colors.push(idx as u32);
        }

        // Build edge_keys + canonical (min, max) pairs.
        let m = source_idx.len();
        let mut edge_keys: AHashSet<u64> = AHashSet::with_capacity(m * 2);
        let mut canon: Vec<(u32, u32)> = Vec::with_capacity(m);
        for k in 0..m {
            let s = source_idx[k];
            let t = target_idx[k];
            if s == t {
                continue; // self-loop, drop
            }
            let (lo, hi) = if s < t { (s, t) } else { (t, s) };
            if (lo as usize) >= n || (hi as usize) >= n {
                return Err(format!(
                    "edge endpoint out of range: ({}, {}); n_nodes = {}",
                    lo, hi, n
                ));
            }
            let key = pack_edge(lo, hi);
            if edge_keys.insert(key) {
                canon.push((lo, hi));
            }
        }

        // Build CSR: count degrees, then fill.
        let mut degrees: Vec<u32> = vec![0; n];
        for &(lo, hi) in &canon {
            degrees[lo as usize] += 1;
            degrees[hi as usize] += 1;
        }
        let mut neigh_offsets: Vec<u32> = Vec::with_capacity(n + 1);
        neigh_offsets.push(0);
        for d in &degrees {
            let last = *neigh_offsets.last().unwrap();
            neigh_offsets.push(last + d);
        }
        let total = *neigh_offsets.last().unwrap() as usize;
        let mut neighbors: Vec<u32> = vec![u32::MAX; total];
        let mut cursors: Vec<u32> = neigh_offsets[..n].to_vec();
        for &(lo, hi) in &canon {
            let cu = cursors[lo as usize] as usize;
            neighbors[cu] = hi;
            cursors[lo as usize] += 1;
            let cv = cursors[hi as usize] as usize;
            neighbors[cv] = lo;
            cursors[hi as usize] += 1;
        }

        // Sort each vertex's neighbor slice. Required by the size-3 closed
        // (sorted-merge triangle) and size-4 algorithms.
        for v in 0..n {
            let lo = neigh_offsets[v] as usize;
            let hi = neigh_offsets[v + 1] as usize;
            neighbors[lo..hi].sort_unstable();
        }

        Ok(SpatialGraphRs {
            n_nodes: n,
            cell_ids,
            neigh_offsets,
            neighbors,
            colors,
            color_levels: levels,
            edge_keys,
        })
    }

    #[inline]
    pub fn neighbors_of(&self, v: u32) -> &[u32] {
        let lo = self.neigh_offsets[v as usize] as usize;
        let hi = self.neigh_offsets[v as usize + 1] as usize;
        &self.neighbors[lo..hi]
    }

    #[inline]
    pub fn has_edge(&self, a: u32, b: u32) -> bool {
        self.edge_keys.contains(&pack_edge(a, b))
    }

    #[inline]
    pub fn color_of(&self, v: u32) -> &str {
        &self.color_levels[self.colors[v as usize] as usize]
    }
}
