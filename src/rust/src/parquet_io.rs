//! Direct parquet → `SpatialGraphRs` builder.
//!
//! Reads `nodes.parquet` (cell_id, x, y, z, cell_type, sample_id,
//! region_id, …) and `edges.parquet` (source, target, distance,
//! sample_id) via `arrow-rs` + `parquet` crates and constructs the
//! in-memory CSR adjacency without round-tripping the bulk data
//! through R. Only string columns we need (`cell_id`, `cell_type`,
//! `source`, `target`) are materialized.

use std::fs::File;

use ahash::AHashMap;
use arrow::array::{Array, StringArray};
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
