# 140k GFP Structure Prediction and 3Di Pipeline

This directory contains three production-oriented modules:

- `generate_a3m.py`: replace the WT A3M query sequence with each mutant sequence and split outputs into `batch_0` ... `batch_4`.
- `run_colabfold_cluster.sh`: launch five `colabfold_batch` workers bound to `CUDA_VISIBLE_DEVICES=0..4`.
- `extract_3di.py`: build a Foldseek database from ColabFold PDB outputs and export amino-acid + 3Di sequences.

## 1. Generate Mutant A3M Files

```bash
python pipeline/generate_a3m.py \
  --wt-a3m wt.a3m \
  --mutants-fasta mutants.fasta \
  --out-dir a3m_batches \
  --num-batches 5 \
  --skip-existing
```

By default, only the first/query sequence in `wt.a3m` is replaced; the original query header is preserved. If your ColabFold version uses the A3M query header for output naming, add:

```bash
--rename-query-header
```

## 2. Run LocalColabFold on 5 GPUs

```bash
bash pipeline/run_colabfold_cluster.sh a3m_batches colabfold_outputs logs_colabfold
```

Fixed ColabFold parameters:

- `--msa-mode custom`
- `--num-recycle 3`
- `--num-relax 0`
- `--num-models 1`
- `--model-order 1`
- `--stop-at-score 100`
- `XLA_PYTHON_CLIENT_MEM_FRACTION=0.9`

If `colabfold_batch` is not in `PATH`:

```bash
COLABFOLD_BIN=/path/to/colabfold_batch bash pipeline/run_colabfold_cluster.sh a3m_batches colabfold_outputs logs_colabfold
```

## 3. Extract Foldseek 3Di

```bash
python pipeline/extract_3di.py \
  --pdb-dir colabfold_outputs \
  --out-csv output_3di.csv \
  --work-dir foldseek_work \
  --threads 32
```

Output CSV columns:

- `Sequence_ID`
- `Amino_Acid_Seq`
- `3Di_Seq`
