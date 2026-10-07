//! Direct parquet → `SpatialGraphRs` builder.
//!
//! Reads `nodes.parquet` (cell_id, x, y, z, cell_type, sample_id,
//! region_id, …) and `edges.parquet` (source, target, distance,
//! sample_id) via `arrow-rs` + `parquet` crates and constructs the
//! in-memory CSR adjacency without round-tripping the bulk data
//! through R. Only string columns we need (`cell_id`, `cell_type`,
//! `source`, `target`) are materialized.

use std::fs::File;

use ahash::{AHashMap, AHashSet};
use arrow::ffi_stream::{ArrowArrayStreamReader, FFI_ArrowArrayStream};
use arrow::array::{Array, Int32Array, Int64Array, StringArray};
use arrow::datatypes::DataType;
use arrow::record_batch::RecordBatch;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

use crate::graph::SpatialGraphRs;

fn read_parquet(path: &str) -> Result<Vec<RecordBatch>, String> {
    let file = File::open(path).map_err(|e| format!("open {}: {}", path, e))?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file)
        .map_err(|e| format!("parquet builder for {}: {}", path, e))?;
    let reader = builder.build()
        .map_err(|e| format!("parquet reader for {}: {}", path, e))?;
    let mut out = Vec::new();
    for batch in reader {
        out.push(batch.map_err(|e| format!("parquet batch read: {}", e))?);
    }
    Ok(out)
}

fn col_string<'a>(batch: &'a RecordBatch, name: &str) -> Result<&'a StringArray, String> {
    let col = batch.column_by_name(name)
        .ok_or_else(|| format!("missing column '{}'", name))?;
    col.as_any().downcast_ref::<StringArray>()
        .ok_or_else(|| format!("column '{}' is not a string array", name))
}

pub fn graph_from_parquet(
    nodes_path: &str,
    edges_path: &str,
) -> Result<SpatialGraphRs, String> {
    let node_batches = read_parquet(nodes_path)?;
    let edge_batches = read_parquet(edges_path)?;

    let mut cell_ids: Vec<String> = Vec::new();
    let mut cell_types: Vec<String> = Vec::new();
    for b in &node_batches {
        let id_arr = col_string(b, "cell_id")?;
        let ct_arr = col_string(b, "cell_type")?;
        for i in 0..b.num_rows() {
            cell_ids.push(id_arr.value(i).to_string());
            cell_types.push(ct_arr.value(i).to_string());
        }
    }

    // Build cell_id → index map for edge resolution.
    let mut id_to_idx: AHashMap<String, u32> = AHashMap::with_capacity(cell_ids.len());
    for (i, id) in cell_ids.iter().enumerate() {
        id_to_idx.insert(id.clone(), i as u32);
    }

    let mut sources: Vec<u32> = Vec::new();
    let mut targets: Vec<u32> = Vec::new();
    for b in &edge_batches {
        let src = col_string(b, "source")?;
        let tgt = col_string(b, "target")?;
        for i in 0..b.num_rows() {
            let s = src.value(i);
            let t = tgt.value(i);
            let si = *id_to_idx.get(s)
                .ok_or_else(|| format!("edge source '{}' not in nodes", s))?;
            let ti = *id_to_idx.get(t)
                .ok_or_else(|| format!("edge target '{}' not in nodes", t))?;
            sources.push(si);
            targets.push(ti);
        }
    }

    SpatialGraphRs::build(cell_ids, cell_types, &sources, &targets)
}


// --- GiottoDisk parquetEdgeStore -----------------------------------------
//
// GiottoDisk stores networks with the node ids already interned:
//
//   edges/  from_id int32|64, to_id int32|64, weight float32, distance float32
//   nodes/  row_index, node_id string, int_id int32|64
//
// So the string -> index hashing that `graph_from_parquet` has to do is work
// GiottoDisk already did at write time. Reading the int columns straight
// through skips it entirely: no AHashMap, no String per cell.
//
// Cell type labels deliberately do not live in the edge store -- the node
// sidecar carries ids only -- so they are supplied by the caller, aligned to
// `int_id`.

/// Read an integer column as `u32`, accepting either int32 or int64 because
/// GiottoDisk promotes to int64 past 2^31 nodes.
fn col_ints(batch: &RecordBatch, name: &str) -> Result<Vec<u32>, String> {
    let col = batch
        .column_by_name(name)
        .ok_or_else(|| format!("missing column '{}'", name))?;
    match col.data_type() {
        DataType::Int32 => {
            let a = col
                .as_any()
                .downcast_ref::<Int32Array>()
                .ok_or_else(|| format!("column '{}' is not int32", name))?;
            Ok((0..a.len()).map(|i| a.value(i) as u32).collect())
        }
        DataType::Int64 => {
            let a = col
                .as_any()
                .downcast_ref::<Int64Array>()
                .ok_or_else(|| format!("column '{}' is not int64", name))?;
            Ok((0..a.len()).map(|i| a.value(i) as u32).collect())
        }
        other => Err(format!(
            "column '{}' has type {:?}; expected int32 or int64",
            name, other
        )),
    }
}

/// Node ids and their integer codes from a `parquetEdgeStore` node sidecar.
pub fn nodes_from_edge_store(nodes_path: &str) -> Result<(Vec<String>, Vec<u32>), String> {
    let batches = read_parquet(nodes_path)?;
    let mut ids = Vec::new();
    let mut ints = Vec::new();
    for b in &batches {
        let id = col_string(b, "node_id")?;
        let iv = col_ints(b, "int_id")?;
        for i in 0..b.num_rows() {
            ids.push(id.value(i).to_string());
            ints.push(iv[i]);
        }
    }
    Ok((ids, ints))
}

/// Build a graph directly from a GiottoDisk `parquetEdgeStore`.
///
/// `colors` must be one code per node, ordered by the sidecar's `int_id`.
/// Edge endpoints are remapped from `int_id` space onto a dense `0..n` range,
/// so an id universe with gaps (a subset store) works unchanged.
pub fn graph_from_edge_store(
    nodes_path: &str,
    edges_path: &str,
    colors: Vec<u32>,
    color_levels: Vec<String>,
) -> Result<SpatialGraphRs, String> {
    let (_ids, int_ids) = nodes_from_edge_store(nodes_path)?;
    let n = int_ids.len();
    if colors.len() != n {
        return Err(format!(
            "colors length {} != node count {} in {}",
            colors.len(),
            n,
            nodes_path
        ));
    }
    // int_id -> dense position; the sidecar is not required to be 0-based or
    // contiguous, and a subset store is not
    let mut dense: AHashMap<u32, u32> = AHashMap::with_capacity(n);
    for (pos, &iid) in int_ids.iter().enumerate() {
        dense.insert(iid, pos as u32);
    }

    let edge_batches = read_parquet(edges_path)?;
    let mut sources = Vec::new();
    let mut targets = Vec::new();
    for b in &edge_batches {
        let f = col_ints(b, "from_id")?;
        let t = col_ints(b, "to_id")?;
        for i in 0..b.num_rows() {
            let (a, c) = (f[i], t[i]);
            let si = *dense
                .get(&a)
                .ok_or_else(|| format!("edge from_id {} not in the node sidecar", a))?;
            let ti = *dense
                .get(&c)
                .ok_or_else(|| format!("edge to_id {} not in the node sidecar", c))?;
            sources.push(si);
            targets.push(ti);
        }
    }
    SpatialGraphRs::build_coded(n, colors, color_levels, &sources, &targets)
}

/// Build a graph from an Arrow stream of edges plus the node sidecar R holds.
///
/// The counterpart of [`graph_from_edge_store`] for a producer that is not a
/// file: R applies whatever narrowing it owes -- a `parquetEdgeStore`'s
/// pending `@ops`, a partitioned dataset, a query over something that is not
/// parquet at all -- and hands over the result as an Arrow stream. Reading
/// the store's files by path cannot see any of that, because a pending subset
/// is not on disk.
///
/// `stream_addr` is the address of a populated `ArrowArrayStream` struct whose
/// ownership moves here; the reader releases it on drop. The stream must yield
/// integer `from_id` / `to_id` columns.
///
/// The node set is the stream's endpoints -- the producer decides which cells
/// are in the network. `int_ids` / `colors` are a label lookup, not the node
/// set: they may cover more ids than the stream uses, and an endpoint missing
/// from them is an error. Nodes are numbered in ascending `int_id` order, not
/// arrival order, because batch order from a multi-file dataset is not fixed
/// and the null permutes over node positions.
pub fn graph_from_edge_stream(
    stream_addr: usize,
    int_ids: Vec<u32>,
    colors: Vec<u32>,
    color_levels: Vec<String>,
) -> Result<SpatialGraphRs, String> {
    if colors.len() != int_ids.len() {
        return Err(format!(
            "colors length {} != int_ids length {}",
            colors.len(),
            int_ids.len()
        ));
    }
    if stream_addr == 0 {
        return Err("the ArrowArrayStream pointer is null".to_string());
    }

    let mut label: AHashMap<u32, u32> = AHashMap::with_capacity(int_ids.len());
    for (&iid, &c) in int_ids.iter().zip(colors.iter()) {
        if label.insert(iid, c).is_some() {
            return Err(format!("int_id {} appears more than once in the label lookup", iid));
        }
    }

    let raw = stream_addr as *mut FFI_ArrowArrayStream;
    let reader = unsafe { ArrowArrayStreamReader::from_raw(raw) }
        .map_err(|e| format!("could not import the Arrow stream: {}", e))?;

    // Endpoints are collected as raw ids and renumbered in place once the
    // node set is known, so peak memory stays at one u32 per endpoint plus a
    // per-node map -- the same as translating on arrival.
    let mut sources: Vec<u32> = Vec::new();
    let mut targets: Vec<u32> = Vec::new();
    let mut seen: AHashSet<u32> = AHashSet::new();
    for batch in reader {
        let b = batch.map_err(|e| format!("could not pull an Arrow batch: {}", e))?;
        let f = col_ints(&b, "from_id")?;
        let t = col_ints(&b, "to_id")?;
        seen.extend(f.iter().copied());
        seen.extend(t.iter().copied());
        sources.extend(f);
        targets.extend(t);
    }

    let mut nodes: Vec<u32> = seen.into_iter().collect();
    nodes.sort_unstable();
    let n = nodes.len();
    let mut dense: AHashMap<u32, u32> = AHashMap::with_capacity(n);
    let mut node_colors: Vec<u32> = Vec::with_capacity(n);
    for (pos, &iid) in nodes.iter().enumerate() {
        let c = *label
            .get(&iid)
            .ok_or_else(|| format!("edge endpoint int_id {} has no label", iid))?;
        dense.insert(iid, pos as u32);
        node_colors.push(c);
    }
    drop(nodes);
    drop(label);
    for v in sources.iter_mut().chain(targets.iter_mut()) {
        *v = dense[v];
    }
    SpatialGraphRs::build_coded(n, node_colors, color_levels, &sources, &targets)
}
