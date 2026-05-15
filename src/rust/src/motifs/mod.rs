//! Motif-result data structures shared across size-2/3/4 kernels.
//!
//! These mirror the R-side `Motifs` schema bit-for-bit:
//!   * `catalog`      → data.table(motif_id, size, canonical_iso_class, color_tuple, count)
//!   * `incidence`    → data.table(motif_id, instance_id, cell_id, orbit_id)
//!   * `instance_meta` (per topology class) → wide tables of cell ids in
//!     canonical role order, plus per-position degree for size-4.
//!
//! Each kernel appends rows to one or more of these tables; the final
//! `MotifResult` is shipped back to R wholesale and wrapped as an S3
//! `Motifs` object.

use std::collections::HashMap;

pub mod size2;
pub mod size3;
pub mod size4;

#[derive(Default, Debug)]
pub struct CatalogTable {
    pub motif_id: Vec<String>,
    pub size: Vec<i32>,
    pub canonical_iso_class: Vec<String>,
    pub color_tuple: Vec<Option<String>>,
    pub count: Vec<i32>,
}

#[derive(Default, Debug)]
pub struct IncidenceTable {
    pub motif_id: Vec<String>,
    pub instance_id: Vec<String>,
    pub cell_id: Vec<String>,
    pub orbit_id: Vec<i32>,
}

#[derive(Default, Debug)]
pub struct EdgeMeta {
    pub instance_id: Vec<String>,
    pub v1: Vec<String>,
    pub v2: Vec<String>,
}

#[derive(Default, Debug)]
pub struct TriangleMeta {
    pub instance_id: Vec<String>,
    pub v1: Vec<String>,
    pub v2: Vec<String>,
    pub v3: Vec<String>,
}

#[derive(Default, Debug)]
pub struct WedgeMeta {
    pub instance_id: Vec<String>,
    pub end1: Vec<String>,
    pub end2: Vec<String>,
    pub center: Vec<String>,
}

#[derive(Default, Debug)]
pub struct Size4Meta {
    pub instance_id: Vec<String>,
    pub w1: Vec<String>,
    pub w2: Vec<String>,
    pub w3: Vec<String>,
    pub w4: Vec<String>,
    pub d1: Vec<i32>,
    pub d2: Vec<i32>,
    pub d3: Vec<i32>,
    pub d4: Vec<i32>,
}

#[derive(Default, Debug)]
pub struct MotifResult {
    pub catalog: CatalogTable,
    pub incidence: IncidenceTable,
    pub edge_meta: Option<EdgeMeta>,
    pub triangle_meta: Option<TriangleMeta>,
    pub wedge_meta: Option<WedgeMeta>,
    /// Class name (`"claw"`, `"path"`, `"cycle"`, `"paw"`, `"diamond"`,
    /// `"K4"`) → wide instance table. Populated only by the size-4 kernel.
    pub size4_meta: HashMap<String, Size4Meta>,
}

/// Convert a `(motif_id → count)` map plus a parallel
/// `(motif_id → (size, canonical_iso_class, color_tuple))` lookup into
/// a `CatalogTable`. The lookup is needed because we don't keep size /
/// class info in the count map (they're constant per kernel).
pub fn catalog_from_counts(
    counts: &HashMap<String, i32>,
    info: &HashMap<String, (i32, String, Option<String>)>,
) -> CatalogTable {
    let mut keys: Vec<&String> = counts.keys().collect();
    keys.sort();
    let mut out = CatalogTable::default();
    for k in keys {
        let (sz, cls, ctup) = &info[k];
        out.motif_id.push(k.clone());
        out.size.push(*sz);
        out.canonical_iso_class.push(cls.clone());
        out.color_tuple.push(ctup.clone());
        out.count.push(counts[k]);
    }
    out
}
