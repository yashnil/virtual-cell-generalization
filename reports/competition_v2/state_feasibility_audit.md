# STATE (Arc Institute) feasibility audit for VCC 2026

**Track:** competition track only (Arc Virtual Cell Challenge 2026). Nothing was trained, submitted, or downloaded beyond config/metadata files. The largest file fetched in full was a 16.7 MB `pert_onehot_map.pt`. Parameter counts come from HTTP range reads of checkpoint and safetensors headers.

**Separation from the research track:** this audit does not bear on the frozen research finding that STRING, DepMap and pathway annotation priors failed to externally predict unseen perturbation effects. STATE tests a different hypothesis: that a model pretrained on large amounts of *direct perturbational* single-cell data transfers. A STATE success or failure here must not be read as evidence for or against the annotation-prior result.

**Audit date:** 2026-09-25. Validation contexts A/B/C are live. Final contexts D/E/F are released 2026-10-22 and final submissions are due 2026-11-05.

## Pinned sources

| Item | Identifier |
|---|---|
| github.com/ArcInstitute/state | commit `9bbfe78a434a55205e4de834e1ea99f85f7a3add` (2026-07-23, "raw-count-normalization-h5ad", pyproject version 0.11.3) |
| PyPI `arc-state` | latest 0.11.1 (uploaded 2026-06-29). The repo HEAD (0.11.3) is newer than PyPI. |
| github.com/ArcInstitute/cell-load | commit `9ba45e59f6f8117bb7a21371ad38d67175586d53` |
| HF `arcinstitute/SE-600M` | `5a9a80f44f7ce32ce57059933ef0d735d7c10ce5` (modified 2026-02-24) |
| HF `arcinstitute/SE-100M` | `bc72702320639217128df42673e94a9658f67d24` |
| HF `arcinstitute/ST-HVG-Replogle` | `bb6a9562cbbf1fd152df14cc53b4cc7517c77175` |
| HF `arcinstitute/ST-SE-Replogle` | `e324967ff4cea5ec199e29bcbb5c1f00e5b9d69c` |
| HF `arcinstitute/st-se-replogle-full` | `d1441f5587ace12c46247b3703c7d25f24f8abb6` |
| HF `arcinstitute/st-x-replogle-full` | `48ad5f70215ab4c58caa5a68e77d837601d29d35` |
| HF `arcinstitute/ST-HVG-Tahoe` | `ca6b751972493f8448e3256d1340ae70ad43e1e7` |
| HF `arcinstitute/ST-SE-Tahoe` | `03b1971d7cc93a7535fd2e957c6948dba267378b` |
| HF `arcinstitute/ST-HVG-Parse` | `827df657cd87063026d2abefd73e99036755ac1b` |
| VCC 2025 Colab ("Train STATE for the Virtual Cell Challenge") | Drive id `1QKOtYP7bMpdgDJEipDxaJqOchv7oQ-_l` (fetched 2026-09-25) |
| VCC 2025 support set | `https://storage.googleapis.com/vcc_data_prod/datasets/state/competition_support_set.zip` (8.72 GB, last-modified 2025-11-16; only the zip listing was read) |

Local scratch copies (not in the repo) are under `/private/tmp/claude-502/.../scratchpad/{state,cell-load,hf}`.

## Summary verdict

| Question | Verdict | Confidence |
|---|---|---|
| 1. Relevant pretrained models | ST checkpoints exist only for Replogle-Nadig (4 essential-screen lines, CRISPRi), Tahoe-100M (drugs) and Parse-PBMC (cytokines). SE-600M and SE-100M are embedding models. There is **no public ST checkpoint trained on the VCC 2025 H1 data, and none uses gene-embedding perturbation features.** | Verified (HF API + configs) |
| 2. Checkpoint sizes / params | ST 49–109 M params (0.47–1.1 GB `.ckpt` including optimizer state). SE-600M has 715 M params (2.86 GB safetensors). SE-100M has 217 M params (0.87 GB). | Verified (header range reads) |
| 3. License | Code is CC BY-NC-SA 4.0. Weights and outputs fall under the Arc State Model **Non-Commercial** License. The VCC FAQ explicitly allows State code for any entrant and pretrained checkpoints for **non-commercial entrants**; commercial entrants must request a license (90-day free trial). | Verified (FAQ text) |
| 4. Gene space | Public ST checkpoints output 2,000 HVGs (Replogle/Tahoe/Parse) or 6,546 genes (Replogle "full"). **None outputs all 18,533 genes.** Of the Replogle HVGs, 1,877/2,000 are among the VCC genes, and 6,197/6,546 of the "full" genes. Output is log1p(CP10k), not counts. | Verified |
| 5. Perturbation encoding / unseen targets | Every public ST checkpoint uses `pert_rep: onehot`, `perturbation_features_file: null`. **0/300 VCC 2026 targets are in the Replogle ST one-hot map (2,024 perts).** `state tx infer` silently maps unknown perturbations to the **control vector**, so the pretrained ST cannot represent any VCC target and would output ≈ control, i.e. the context-mean baseline (scaled score ≈ 0). The code *does* support ESM2 gene-embedding features: ESM2 (5,120-d) embeddings exist for **300/300** targets (SE-600M `protein_embeddings.pt`). Using them requires **training a new ST**. | Verified |
| 6. Unseen cell context | Yes by design: ST conditions on a set of control cells, so A/B/C/D/E/F controls go in directly. Caveats: batch-encoder fallback and the 10x Flex vs 10x 3′ platform shift. | Verified (code); transfer quality unverified |
| 7. Compute | Inference of a ~50 M-param ST over 360k cells is trivial (minutes; CPU feasible on M1). SE-600M embedding of 55k controls on M1 CPU is slow (hours, estimate). **Training a new ESM2-featurized ST (the only route that can represent the targets) needs a CUDA GPU.** The 2025 Colab ran at about 1.2 it/s on a T4, so 40k steps is about 9 T4-hours. | Inference: verified path; timings estimated |
| 8. Counts output | No. Outputs are log1p-normalized expression, clipped to [0, 14]. A back-transform (expm1, library-size rescale, stochastic integerisation) plus gap-filling of non-modelled genes is required. | Verified |
| 9. VCC 2026 rules | Any model and any data are allowed if you have the rights. Finalists must disclose published models used. State/Stack use is explicitly addressed (row 3). | Verified (virtualcellchallenge.org JS bundle, arcinstitute.org news) |

**Overall:** the **off-the-shelf pretrained STATE is not feasible** for VCC 2026. Its checkpoints cannot represent any of the 300 targets, and inference would reduce to a control/context-mean prediction. STATE is only usable as a *framework*: train a new ST with ESM2 perturbation features on public CRISPRi Perturb-seq, conditioned on VCC controls. That is a GPU training project of roughly 10–30 GPU-hours, not an inference-only test, and it falls outside this audit's no-training remit.

---

## 1. Pretrained models and training data

HF organisation listing (`curl https://huggingface.co/api/models?author=arcinstitute`) shows these STATE-family models: SE-600M, SE-100M, ST-HVG-Replogle, ST-SE-Replogle, st-se-replogle-full, st-x-replogle-full, ST-HVG-Tahoe, ST-SE-Tahoe, ST-HVG-Parse and ST-SE-Parse. The same listing includes Stack-Large and Stack-Large-Aligned, which are a separate Arc model.

| Model | Trained on (from config `toml_config_path` / pert map) | Perturbation vocabulary | Splits |
|---|---|---|---|
| ST-HVG-Replogle / ST-SE-Replogle | Replogle-Nadig essential CRISPRi screens: K562, RPE1, HepG2, Jurkat (`/large_storage/.../replogle_nogwps_v2/hepg2_zeroshot.toml`. The "nogwps" path suggests the K562 genome-wide screen was excluded; inferred from path name, unverified.) | 2,024 genes, one-hot (`pert_onehot_map.pt`, verified) | `zeroshot/{cell}` = that cell line held out; `fewshot/{cell}` |
| st-se-replogle-full / st-x-replogle-full | Replogle-Nadig, `output_space: all` over 6,546 genes, `downsample: 0.99` | same 2,024 one-hot | `{cell}_0.99` |
| ST-HVG-Tahoe / ST-SE-Tahoe | Tahoe-100M (`/data/tahoe_se/generalization_zeroshot.toml`. Zero-shot test lines: C32, HOP62, HepG2/C3A, Hs 766T, PANC-1.) | 1,138 drug×dose one-hot (e.g. `[('(S)-Crizotinib', 0.5, 'uM')]`) | zeroshot / fewshot |
| ST-HVG-Parse / ST-SE-Parse | Parse PBMC cytokines | 91 cytokines one-hot (e.g. `4-1BBL`, `BAFF`) | 5 splits × zeroshot/fewshot. HF README: "Certain uses of a model trained with Parse data may require a license from Parse Biosciences, Inc." |
| SE-600M / SE-100M | Observational scRNA-seq. Arc news: SE trained on "167 million" human cells; ST on ">100 million perturbed cells across 70 cell contexts" (https://arcinstitute.org/news/virtual-cell-model-state, via search snippet; paper PDF fetch was rate-limited, so treat as secondary-source). | n/a (gene tokens = ESM2 protein embeddings, 19,790 genes × 5,120-d) | n/a |

**No public ST checkpoint was trained on VCC 2025 H1 data.** The 2025 Colab *trains* one from scratch; it does not ship one.

## 2. Checkpoint availability, sizes and parameter counts

All repos are public and **not gated** (`gated: False`). Parameter counts come from `state_dict` shapes (range-read `data.pkl` inside the Lightning ckpt) and safetensors headers.

| Checkpoint | File size | Params (exact) | Notes |
|---|---|---|---|
| ST-HVG-Replogle `zeroshot/hepg2/checkpoints/final.ckpt` | 471.7 MB | 49,396,728 | hidden 328, 8-layer bidirectional Llama (42.7 M), cell_set_len 64, batch_encoder (56 gem groups), output 2,000 HVG |
| ST-SE-Replogle `.../final.ckpt` | 510.1 MB | ≈ 50–53 M (not range-read) | input `X_state` (2,058-d SE embedding), output 2,000 HVG |
| st-x-replogle-full `hepg2_0.99/.../final.ckpt` | 664.2 MB | 70,087,272 | input/output 6,546 genes (`output_space: all`) |
| st-se-replogle-full | 550.9 MB | not range-read | SE input, 6,546-gene output |
| ST-HVG-Tahoe `zeroshot/.../final.ckpt` | 1.07 GB | 108,690,592 | hidden 768, cell_set_len 256 |
| ST-SE-Tahoe | 1.11 GB | not range-read | |
| ST-HVG-Parse / ST-SE-Parse | 540 MB / 579 MB | not range-read | hidden 384, cell_set_len 512 |
| SE-600M `model.safetensors` | 2.86 GB | 715,062,442 (F32), of which 101.3 M is the ESM2 table | also `se600m_epoch16.ckpt` (11.5 GB), `se600m_epoch4.safetensors` (2.82 GB) and `protein_embeddings.pt` (410.9 MB) |
| SE-100M `model.safetensors` | 869 MB | 217,284,776 (F32), incl. 101.3 M ESM2 table | plus `protein_embeddings.pt` (411.1 MB) |

Each ST run directory also contains `config.yaml`, `var_dims.pkl` (gene_names, pert_names), `pert_onehot_map.pt` (16.7 MB for Replogle), `batch_onehot_map.pkl` and `cell_type_onehot_map.pkl`. The HF repos also host large eval artefacts (e.g. `adata_pred.h5ad`, 573 MB), which should be skipped via `allow_patterns`.

## 3. Licensing

Sources are the repo `LICENSE`, `MODEL_LICENSE.md` and `README.md` §Licenses, and the VCC FAQ.

- **Code:** CC BY-NC-SA 4.0 (README: "State code is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0").
- **Weights and outputs:** "The model weights and output are licensed under the Arc Research Institute State Model Non-Commercial License" (README). Key clauses in `MODEL_LICENSE.md` (last updated 2025-06-23):
  - §1.4 defines Non-Commercial Purpose as "not undertaken for direct or indirect monetary compensation…". Research sponsored by a Commercial Entity is *not* non-commercial.
  - §1.5 defines Output as the results of running the model. §2.1 grants creation of Outputs "solely for Non-Commercial Purposes".
  - §3 requires any redistributed Output to carry a citation to Adduri et al. 2025.
  - §4.1 bars use for anything other than Non-Commercial Purposes without a separate license (legal@arcinstitute.org).
- **VCC 2026 FAQ** ("Can I use Arc's State and/or Stack model in the Challenge?"), verbatim from the virtualcellchallenge.org JS bundle:
  > "Any entrant may use Arc's State and/or Stack model code in order to participate in the Challenge (and such use will be considered a Non-Commercial Purpose solely to the extent of your use for participation in the Challenge)… With respect to Arc's State and/or Stack model pretrained checkpoints, you may use them if you are a non-commercial entrant, and if you are a commercial entrant, you may request a commercial license to certain of such pretrained checkpoints, which, if confirmed, includes a 90-day free trial permitting you to use such models in VCC."
- **Verdict:**
  - An individual or academic entrant not acting for a company may use both code and checkpoints, and submitting their outputs is fine.
  - Code-only use (training your own ST) is allowed for **any** entrant.
  - Parse-trained checkpoints carry an extra Parse Biosciences caveat.
  - Finalists must disclose "the model or models – including any published models or publicly available models – you used" (Rules, finalists section).
  - The 2025 Colab asks users to "mention it in your model description".
  - *Not legal advice.* The user's entrant status (personal vs employer) should be confirmed.

## 4. Gene coverage and input space

- `state tx preprocess_train` does `normalize_total` → `log1p` → `highly_variable_genes(n_top_genes)` and stores `.obsm["X_hvg"]` (README; `_preprocess_train.py`).
- ST-HVG-* checkpoints take and emit 2,000 HVG values in **log1p(CP10k)** space. The Replogle HVG list overlaps VCC `gene_names.csv` at **1,877/2,000**, so 123 inputs would be zero-filled. It also covers only about 10% of the 18,533 VCC readout genes.
- ST-SE-* checkpoints take SE-600M embeddings (`X_state`, 2,058-d) and still emit 2,000 HVGs (`output_space: gene`).
- st-x/st-se-replogle-full emit 6,546 genes (6,197 in the VCC gene list).
- SE-600M's gene vocabulary is the ESM2 table of 19,790 genes, which covers 18,190/18,533 VCC genes. SE-100M covers 17,930.
- The 2025 Colab trained an ST with `output_space: all` on 18,080 genes. That is the only demonstrated transcriptome-wide configuration, and it has to be trained.
- Mapping back to 18,533 raw counts requires imputing non-modelled genes, e.g. from the context's control mean or sampled control cells, and inverting the normalisation (§8).

## 5. Perturbation encoding and unseen targets (crucial)

- The encoder (`state_transition.py`) is `pert_encoder = build_mlp(in_dim=pert_dim, …)` applied to `batch["pert_emb"]` and added to the basal (control) cell embeddings before the transformer. `pert_emb` comes from `pert_onehot_map.pt`.
- **Every public ST checkpoint uses one-hot encodings.** All audited configs have `pert_rep: onehot` and `perturbation_features_file: null`. The Replogle map contains 2,024 one-hot vectors (binary, sum 1).
- **Coverage check against VCC 2026 (`data/raw/arc2026/controls/pert_counts.csv`, 300 targets):**

  | Vocabulary | VCC-300 overlap |
  |---|---|
  | ST-Replogle one-hot map (2,024) | **0 / 300** |
  | Replogle K562 essential (scPerturBench copy, 1,972) | 0 |
  | Replogle RPE1 essential (2,017) | 0 |
  | Nadig HepG2 (1,819) | 0 |
  | Nadig Jurkat (2,138) | 0 |
  | VCC 2025 H1 training perts (150) / arch1 (151) | 13 |
  | Kaden 2025 RPE1 (1,837) | **80** |
  | Wessels 2023 (158) | 0 |
  | ESM2 table in SE-600M `protein_embeddings.pt` (19,790) | **300 / 300** |
  | ESM2 table in SE-100M (19,790, different build) | 299 / 300 (missing C5orf22) |
  | Replogle K562 genome-wide (GWPS, in VCC 2025 support set as `k562_gwps.h5`) | **UNVERIFIED** (likely high, since GWPS targets ~9.9k expressed genes; not read because the file is 8 GB) |

  VCC 2026 appears to have deliberately excluded the Replogle-Nadig essential-gene panel.
- **Behaviour on unknown perturbations (`_infer.py` L798–806):** `vec = pert_onehot_map.get(map_key, None)`; if None it uses `default_pert_vec`, which is the **control one-hot**, and prints only "pert '…' not in mapping; using control fallback one-hot." So running a pretrained ST on VCC targets **silently returns control-conditioned output for every target**. That is effectively the context-mean baseline, and it scores ≈ 0 by construction of the VCC scaled metric.
- **Zero-shot unseen perturbation via gene embeddings:** cell-load supports `perturbation_features_file` (`perturbation_dataloader.py` L567–585). It loads a `{pert_name: tensor}` dict, zero-fills any missing pert, and saves the dict as the run's `pert_onehot_map.pt`. The VCC 2025 Colab uses exactly this with `competition_support_set/ESM2_pert_features.pt` (410,886,729 bytes). That size is byte-identical to SE-600M `protein_embeddings.pt`, so it is very likely the same 19,790×5,120 ESM2 table (inferred, not hash-checked). The Colab comment reads: "using ESM2 featurizations of genes as the perturbation embeddings. Note that we are now generalizing across both contexts and perturbations."
- **Conclusion:** zero-shot representation of all 300 targets is possible **only with a newly trained ESM2-featurized ST**. It is impossible with any shipped checkpoint, because pert_dim is 2,024, 1,138 or 91 one-hot and cannot accept 5,120-d features.

  The Colab itself cautions: "STATE was not designed specifically for this setting (1. transcriptome-wide effect prediction and 2. predicting effects of perturbations unseen in any context in the test dataset)". ESM2 protein embeddings are a sequence-similarity prior. Whether they carry transferable knockdown-response signal is exactly the kind of prior our research track found weak for STRING/DepMap/pathways. That is a separate hypothesis, but it lowers expectations.

## 6. Unseen cell context

ST is a set-to-set model. It samples a set of control cells from the target context (`cell_set_len` = 64/256/512), adds the perturbation embedding, and predicts the perturbed set. `state tx infer --all-perts` or `--tsv` clones control cells as "virtual" perturbed cells, so only controls are needed. This matches VCC's inputs.

Caveats:

- **Batch encoder:** Replogle checkpoints have `batch_encoder: True` (56 gem groups). VCC contexts have no matching label, so infer falls back to batch index 0 ("using index 0 as fallback") or disables the encoder if no batch column exists.
- **Platform shift:** VCC 2026 uses **10x Flex** (probe-based; FAQ) at a median of about 20k UMIs per cell. Replogle and Nadig used 10x 3′.
- **Context-transfer evidence is limited:** the published Replogle "zeroshot" checkpoints hold out one of only 4 lines, and the perturbation was *seen* in the other lines.

## 7. Compute requirements

| Step | Size | Hardware | Estimate | Basis |
|---|---|---|---|---|
| ST inference, 900 groups × 400 cells = 360k cells | ~50–110 M params, 64–512-cell sets, ≈ 5.6k forward passes at set 64 | M1 CPU | minutes (estimate) | `_infer.py` has no hard CUDA dependency: it uses `next(model.parameters()).device`, and Lightning loads to CPU when CUDA is absent (standard behaviour, not executed here). The Tahoe Colab did 1.84 M cells in 3 min 51 s on a Colab GPU (notebook output). |
| RAM for HVG inference | 120k virtual cells × 2,000 × 4 B ≈ 1 GB per context | M1 16 GB | OK | |
| RAM for an 18k-gene `output_space: all` ST | 120k × 18,533 × 4 B ≈ 8.9 GB dense per context | M1 16 GB | tight; do one context at a time or `.npy` output | `_infer.py` materialises dense `sim_X` |
| SE-600M embedding of 55,200 controls (needed only for ST-SE-*) | 715 M params; code autocasts only for `cuda`, else CPU | M1 CPU | several hours (unverified estimate); T4 roughly 10–30 min | `emb/inference.py` L130: `"cuda" if available else "cpu"` (no MPS path) |
| **Training an ESM2-featurized ST** (the only viable route) | ~50 M params, 40k steps | CUDA GPU | T4 ≈ 9 h (Colab log: 1.19 it/s, "Epoch 103 … 60/213"). A10G/L4 ≈ 3–5 h, A100 ≈ 1.5–3 h (scaled estimates) | 2025 Colab (`gpuType: T4`, `machine_shape: hm`) |

**Cost (approximate; prices not verified today):**

- Colab Pro ($9.99 per 100 compute units; a T4 uses about 2 CU/h, an A100 about 12 CU/h) covers one 40k-step T4 run.
- Cloud A10G/L4 at about $0.8–1.2/h gives about $5 per run.
- A realistic sweep of 3–5 runs × 2 feature sets costs about $30–100 and takes 20–50 GPU-hours.
- Mac M1 training is not realistic: there is no CUDA, `state tx train` configs assume `device: cuda`, and MPS support is untested.

## 8. Counts output

- Outputs are log-normalized. `infer` writes `adata.X` or `obsm['X_hvg']` in log1p(CP10k) space, clipped to [0, 14] (`clip_array`, `_infer.py` L147).
- A `pert_cell_counts_preds` path exists only for `nb_decoder`/count heads. All shipped checkpoints have `nb_decoder: False`.
- The July 2026 commit `13bc762` ("normalize raw counts in model workflows") adds `log1p_from_raw_counts`, which normalizes raw-count *inputs* inside the model. Outputs are still log space.
- The 2025 pipeline submitted log-normalized values (cell-eval prep: "Input is found to be log-normalized already").
- Required conversion for VCC 2026 raw integer counts:
  1. For modelled genes, compute x̂ = expm1(pred).
  2. Rescale to a per-cell library size drawn from the context's control UMI distribution (median ≈ 20k), i.e. counts_mean = x̂ / 1e4 × L.
  3. For non-modelled genes, use control-cell values from a matched sampled control cell.
  4. Integerise stochastically (Poisson or NB with a dispersion fitted on controls), or round.

  This reuses the existing count-space baseline machinery in this repo (frozen "Arc count-space baseline v1").

## 9. VCC 2026 rules on pretrained models and public data

Sources are virtualcellchallenge.org (FAQ/Rules text extracted from the site's JS bundle, since the pages are client-rendered) and https://arcinstitute.org/news/virtual-cell-challenge-2026.

- "Participants may use any modeling strategy and train their models on any data." (Arc news)
- FAQ "Can I use proprietary models or data…?": "Yes, so long as you have the right to do so."
- FAQ "Can I use other data to train or fine tune my model?": "As long as you have the appropriate permissions to use the data, you may use any data to improve the performance of your predictions. This includes proprietary data…"
- State/Stack: see §3.
- Finalists must publicly describe the datasets used, the published or public models used, and the training, adaptation and inference approach.
- **Important for any approach:** FAQ: "All three contexts in a round share the same perturbation panel; **the validation and final test rounds use different panels**." The D/E/F targets are unknown until 2026-10-22, so any target-specific coverage check must be repeated then.
- The Cell paper (https://www.cell.com/cell/fulltext/S0092-8674(26)00931-1) returned HTTP 403. Its details on target selection are **unverified**.

## Next step (only if pursued; requires lifting the no-training constraint)

An *inference-only* test of a pretrained checkpoint is **not informative**: it would reproduce the control/context-mean baseline because 0/300 targets are in the vocabulary. The cheapest meaningful check is to confirm that fallback in 10 minutes on the Mac, then decide whether to fund a GPU training run.

**A. Local smoke test (M1, CPU, ~0.5 GB download).** This confirms that install works on macOS and that VCC targets hit the control fallback.

```bash
cd /private/tmp/claude-502/-Users-yashnilmohanty-Desktop-virtual-cell-generalization/5d8cf2fe-448e-4e0f-b895-15f358de2559/scratchpad
uv tool install arc-state            # or: cd state && uv run state --help  (repo HEAD 0.11.3)
uv run --with huggingface_hub python - <<'EOF'
from huggingface_hub import snapshot_download
snapshot_download("arcinstitute/ST-HVG-Replogle",
    revision="bb6a9562cbbf1fd152df14cc53b4cc7517c77175",
    allow_patterns=["zeroshot/hepg2/*.pkl","zeroshot/hepg2/*.pt","zeroshot/hepg2/*.yaml",
                    "zeroshot/hepg2/*.torch","zeroshot/hepg2/checkpoints/final.ckpt"],
    local_dir="ST-HVG-Replogle")
EOF
# Build input: context_A controls -> normalize_total+log1p -> reindex to var_dims['gene_names'] (2,000, zero-fill 123)
#   -> obsm['X_hvg']; obs['gene']='non-targeting'; obs['cell_type']='A'; subsample 2,000 controls.
# Then pad 5 VCC targets (TSV columns: perturbation,num_cells) and run:
state tx infer --model-dir ST-HVG-Replogle/zeroshot/hepg2 \
  --checkpoint ST-HVG-Replogle/zeroshot/hepg2/checkpoints/final.ckpt \
  --adata ctxA_hvg.h5ad --pert-col gene --embed-key X_hvg \
  --tsv five_targets.tsv --output ctxA_pred.h5ad
# Expected: "pert 'ACLY' not in mapping; using control fallback one-hot." for all 5.
```

**B. Meaningful experiment (Colab or cloud GPU, 1 run ≈ 9 T4-h or ≈ 3–5 A10G-h, ≈ $5–15).**

1. Build a training TOML from public CRISPRi Perturb-seq whose targets overlap VCC:
   - K562 GWPS: first verify coverage of the 300.
   - Kaden 2025 RPE1: 80/300 overlap.
   - Replogle/Nadig essential lines, to learn the control→perturbed map.
   - VCC 2025 H1: 13 overlap.
2. Hold out VCC-300 targets that appear in the public data as a pseudo-unseen validation set.
3. Train:

   ```
   state tx train model=state data.kwargs.perturbation_features_file=<ESM2 table> data.kwargs.output_space=all data.kwargs.control_pert=non-targeting data.kwargs.pert_col=target_gene training.max_steps=40000
   ```

   Use the gene set = VCC 18,533 ∩ training genes, and `log1p_from_raw_counts=true` if inputs are raw.
4. Evaluate on held-out targets and on held-out *cell lines* before touching the leaderboard.
5. Infer on the A/B/C controls with `--tsv` listing all 300 targets at 400 cells each.
6. Convert to integer counts as in §8, then validate with the cell-eval2 `vcc2026` profile.

Gate: proceed only if held-out-target performance beats the frozen count-space baseline on a local proxy.

## Risks

1. **Target representability (fatal for off-the-shelf use):** 0/300 targets are in any shipped ST vocabulary, and the fallback is silent (control vector), so a naive run looks like it "works".
2. **ESM2 features may not transfer knockdown effects:** the Colab authors state STATE was not designed for unseen-perturbation prediction. This expectation is analogous in kind (though not in hypothesis) to our annotation-prior null.
3. **Final panel differs from the validation panel:** coverage and fine-tuning choices must be re-audited on 2026-10-22 with 14 days left.
4. **Compute and time:** training needs a CUDA GPU (none locally). With about 27 days before the final release and 41 before the deadline, there is room for about 2–4 iteration cycles at most.
5. **Gene-space mismatch:** shipped checkpoints model 2k or 6.5k genes, so the rest must be imputed. A new ST with `output_space: all` at 18.5k genes raises memory needs (about 9 GB per context at inference).
6. **Count conversion:** the log→integer back-transform can distort DE-based metrics (Wilcoxon on counts). This needs local validation.
7. **Platform and batch shift:** 10x Flex vs 3′ data, plus a batch-encoder fallback to index 0.
8. **Licensing:** checkpoints are OK only for a non-commercial entrant. Parse-trained weights carry an extra caveat. Outputs must cite Adduri et al. 2025, and finalists must disclose use.
9. **Version drift:** PyPI 0.11.1 lags repo 0.11.3. Pin the commit SHA `9bbfe78` when installing from source.
10. **Unverified items:** paper training-data numbers (secondary source); ESM2 file identity (size-matched only); GWPS coverage of the VCC-300; M1 wall-clock times; Cell paper details (403).
