# Kailash BP — BigBrainNetworks backbone: read_bigbrain per split, then join_bigbrain
import os

# ---------- Config ----------
outFolder      = config.get("outFolder", "results")
variants_tsv   = config["variants"]                         # TSV: chr, pos, end, variant_id
sumstats_tsv   = config["full_assoc_sumstat"]               # MUST be bgzipped + tabix-indexed
n_splits       = int(config.get("n_splits", 1000))          # Split variants into how many parts?
threads        = int(config.get("threads", 8))
chunk_prefix   = config.get("chunk_prefix", "variants.part")

# MODULES AND SOFTWARE 
R_VERSION = config.get("R_version", "R/4.2.0")
TABIX_VERSION = config.get("TABIX_version", "tabix/0.2.6")

# scripts
split_variants_R   = config.get("split_variants_script", "scripts/split_variants.R")
read_bigbrain_R    = config.get("read_bigbrain_script",  "scripts/read_bigbrain_wrapper.R")
join_bigbrain_R    = config.get("join_bigbrain_script",   "scripts/join_bigbrain_wrapper.R")

# join_bigbrain outputs
combined_inst_gz   = os.path.join(outFolder, "combined_inst_df.tsv.gz")
combined_result_rds= os.path.join(outFolder, "combined_result_df.rds")

# MR settings / outputs
# mr_method   = config.get("mr_method", "ivw")
# mr_min_inst = int(config.get("mr_min_inst", 2))
# mr_dir      = os.path.join(outFolder, "mr")
# os.makedirs(mr_dir, exist_ok=True)

# ---------- Wildcard constraints ----------
wildcard_constraints:
    I = r"\d{5}"

# ---------- Folders ----------
chunks_dir     = os.path.join(outFolder, "chunks")
subset_dir     = os.path.join(outFolder, "subset_sumstats")
rds_dir        = os.path.join(outFolder, "rds")
os.makedirs(chunks_dir, exist_ok=True)
os.makedirs(subset_dir, exist_ok=True)
os.makedirs(rds_dir, exist_ok=True)

# ---------- Derived ----------
pad = lambda i: f"{i:05d}"
chunk_files   = [f"{chunks_dir}/{chunk_prefix}_{i:05d}.bed" for i in range(1, n_splits+1)]
subset_files  = [f"{subset_dir}/{chunk_prefix}_{i:05d}.sumstats.tsv.gz" for i in range(1, n_splits+1)]
rds_files     = [f"{rds_dir}/chunk_{i:05d}.rds" for i in range(1, n_splits+1)]

# ---------- Rules ----------
rule all:
    input:
        combined_inst_gz,
        combined_result_rds

# 1) Split variants into N parts (does not keep header)
rule split_variants:
    input:
        variants_tsv
    output:
        chunk_files
    params:
        script = split_variants_R,   # must write headerless BED chunks
        n      = n_splits,
        outdir = chunks_dir,
        prefix = chunk_prefix
    shell:
        """
        ml {R_VERSION};
        Rscript {params.script} \
          --input {input} \
          --n-splits {params.n} \
          --outdir {params.outdir} \
          --prefix {params.prefix}
        """

# 2) Subset the bgzipped + indexed BigBrain file per split using tabix
rule subset_sumstats:
    input:
        bed      = f"{chunks_dir}/{chunk_prefix}_{{I}}.bed",   # headerless 0-based BED
        sumstats = sumstats_tsv
    output:
        subset   = f"{subset_dir}/{chunk_prefix}_{{I}}.sumstats.tsv.gz"
    params:
        tmp   = lambda wc: f"{subset_dir}/.{chunk_prefix}_{wc.I}.tmp.tsv",      # plain TSV until bgzip
        regs  = lambda wc: f"{subset_dir}/.{chunk_prefix}_{wc.I}.regions",
        tdir  = lambda wc: f"{subset_dir}/.{chunk_prefix}_{wc.I}.parts"
    threads: int(config.get("subset_threads", 8))
    shell:
        r"""
        set -euo pipefail
        ml {TABIX_VERSION}

        # 0) Prep temp dir
        rm -rf {params.tdir}
        mkdir -p {params.tdir}

        # 1) BED(0-based) -> 1-based "chr:start-end"
        awk 'BEGIN{{OFS=""}} {{pos=$2+1; print $1,":",pos,"-",pos}}' {input.bed} > {params.regs}

        # 2) Write header exactly once (respects -S meta-lines from indexing)
        # tabix -H {input.sumstats} > {params.tmp}
        zcat {input.sumstats} | head -n 1 > {params.tmp}

        # 3) If no regions, just compress header and exit
        if [ ! -s {params.regs} ]; then
          bgzip -c {params.tmp} > {output.subset}
          rm -f {params.tmp} {params.regs}
          rmdir {params.tdir}
          exit 0
        fi

        # 4) Split regions into {threads} shards (round-robin) for parallel tabix, no -R
        awk -v T={threads} '{{ fn=sprintf("{params.tdir}/regs.%02d", (NR-1)%T); print > fn }}' {params.regs}

        # 5) For each shard: sequentially query each region -> its own part file (no concurrent appends)
        for f in {params.tdir}/regs.*; do
          [ -s "$f" ] || continue
          of="${{f/regs./part.}}"
          (
            while IFS= read -r r; do
              tabix {input.sumstats} "$r" || true
            done < "$f"
          ) > "$of" &
        done
        wait

        # 6) Concatenate part files in numeric order, then compress
        for p in $(ls {params.tdir}/part.* 2>/dev/null | sort); do
          cat "$p" >> {params.tmp}
        done

        bgzip -c {params.tmp} > {output.subset}

        # 7) Cleanup
        rm -f {params.tmp} {params.regs}
        rm -rf {params.tdir}
        """

# 3) Run read_bigbrain() on each subset; write RDS
rule read_bigbrain_per_split:
    input:
        subset = f"{subset_dir}/{chunk_prefix}_{{I}}.sumstats.tsv.gz"
    output:
        rds    = f"{rds_dir}/chunk_{{I}}.rds"
    params:
        script = read_bigbrain_R
    shell:
        """
        ml {R_VERSION};
        Rscript {params.script} \
          --sumstats {input.subset} \
          --out {output.rds}
        """

# 4) Join all chunk RDS with join_bigbrain(); emit combined outputs
rule join_bigbrain_all:
    input:
        # drive creation of all per-split RDS by expanding the wildcard
        rds = expand(f"{rds_dir}/chunk_{{I}}.rds", I=[f"{i:05d}" for i in range(1, n_splits+1)])
    output:
        inst_gz    = combined_inst_gz,
        result_rds = combined_result_rds
    params:
        script = join_bigbrain_R
    shell:
        """
        ml {R_VERSION};
        Rscript {params.script} \
          --inputs {input.rds} \
          --inst {output.inst_gz} \
          --result {output.result_rds}
        """
