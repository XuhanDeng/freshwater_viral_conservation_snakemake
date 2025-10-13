configfile: "profiles/config.yaml"

rule all:
    input:
        expand("results/10_cluster/1_filtered_seq/{sample}_filterd_final_vircontig.fasta", sample=config["samples"]),
        expand("results/10_cluster/3_clustered/per_sample/{sample}/{sample}_cluster.tsv", sample=config["samples"]),
        expand("results/10_cluster/4_cluster_fasta/{sample}_cluster_representatives.fasta", sample=config["samples"]),
        "results/10_cluster/3_clustered/all_samples_cluster.tsv",
        "results/10_cluster/4_cluster_fasta/all_samples_cluster_representatives.fasta",
        expand("results/11_bowtie2/per_sample/{sample}/{sample}_TPM.tsv", sample=config["samples"]),
        expand("results/11_bowtie2/all_samples/{sample}/{sample}_TPM.tsv", sample=config["samples"])
# Step 1: Quality control with fastp
rule fastp_quality_control:
    input:
        r1="input/{sample}.R1.raw.fastq.gz",
        r2="input/{sample}.R2.raw.fastq.gz"
    output:
        r1_paired="results/1_fastp/{sample}/{sample}_1P.fq.gz",
        r2_paired="results/1_fastp/{sample}/{sample}_2P.fq.gz",
        r1_unpaired="results/1_fastp/{sample}/{sample}_U1.fq.gz",
        r2_unpaired="results/1_fastp/{sample}/{sample}_U2.fq.gz",
        html="results/1_fastp/{sample}/{sample}.fastp.html",
        json="results/1_fastp/{sample}/{sample}.fastp.json"
    conda: "envs/fastp.yaml"
    log:
        out="log/1_fastp/{sample}.log",
        err="log/1_fastp/{sample}.err"
    threads: config["fastp"]["threads"]
    resources:
        slurm_partition=config["regular_partition"],
        runtime = config["runtime"],
        mem_mb_per_cpu = config["regular_memory"],
        cpus_per_task = config["fastp"]["threads"],
        slurm_account = config["slurm_account"]

    params:
        q_cutoff=config["fastp"]["qualified_quality_phred"],
        length_required=config["fastp"]["length_required"]
    shell:
        """
        mkdir -p results/1_fastp/{wildcards.sample}
        mkdir -p log/1_fastp
        fastp --thread {threads} \
            --in1 {input.r1} --in2 {input.r2} \
            --out1 {output.r1_paired} --out2 {output.r2_paired} \
            --unpaired1 {output.r1_unpaired} --unpaired2 {output.r2_unpaired} \
            -h {output.html} -j {output.json} \
            --trim_poly_g --trim_poly_x \
            --qualified_quality_phred {params.q_cutoff} --length_required {params.length_required} \
            --dont_overwrite \
            > {log.out} 2> {log.err}
        """
# Step 2: Build bowtie2 index for T7 reference
rule build_bowtie2_index:
    input:
        ref=config["bowtie2_index"]
    output:
        flag="results/T7_ref/build_index/T7_index.done"
    conda: "envs/bowtie2.yaml"
    log:
        out="log/2_bowtie2_index/build_index.log",
        err="log/2_bowtie2_index/build_index.err"
    threads: config["bowtie2_build"]["threads"]
    resources:
        slurm_partition=config["regular_partition"],
        runtime=config["runtime"],
        mem_mb_per_cpu=config["regular_memory"],
        cpus_per_task=config["bowtie2_build"]["threads"],
        slurm_account=config["slurm_account"]
    shell:
        """
        mkdir -p results/T7_ref/build_index
        mkdir -p log/2_bowtie2_index
        bowtie2-build {input.ref} results/T7_ref/build_index/T7 > {log.out} 2> {log.err}
        touch {output.flag}
        """

# Step 3: Remove T7 phage sequences using bowtie2 and samtools
rule remove_t7_sequences:
    input:
        r1="results/1_fastp/{sample}/{sample}_1P.fq.gz",
        r2="results/1_fastp/{sample}/{sample}_2P.fq.gz",
        idx_flag="results/T7_ref/build_index/T7_index.done"
    output:
        r1_clean="results/2_T7_removal/5_removed_sequence/{sample}/{sample}_host_removed_R1.fastq.gz",
        r2_clean="results/2_T7_removal/5_removed_sequence/{sample}/{sample}_host_removed_R2.fastq.gz"
    conda: "envs/bowtie2_samtools.yaml"
    log:
        out="log/3_T7_removal/{sample}_T7_removal.log",
        err="log/3_T7_removal/{sample}_T7_removal.err"
    threads: config["bowtie2_removal"]["threads"]
    resources:
        slurm_partition=config["regular_partition"],
        runtime=config["runtime"],
        mem_mb_per_cpu=config["regular_memory"],
        cpus_per_task=config["bowtie2_removal"]["threads"],
        slurm_account=config["slurm_account"]
    params:
        sort_mem=lambda wildcards, resources: f"{int(resources.mem_mb_per_cpu * resources.cpus_per_task)}M"
    shell:
        """
        mkdir -p results/2_T7_removal/1_sam/{wildcards.sample} \
                 results/2_T7_removal/2_bam/{wildcards.sample} \
                 results/2_T7_removal/3_unmapped/{wildcards.sample} \
                 results/2_T7_removal/4_unmapped_sorted/{wildcards.sample} \
                 results/2_T7_removal/5_removed_sequence/{wildcards.sample} \
                 log/3_T7_removal
        
        # Initialize log files
        echo "T7 removal started for sample {wildcards.sample}" > {log.out}
        echo "T7 removal started for sample {wildcards.sample}" > {log.err}
        
        # Align reads to T7 reference
        echo "Step 1: Aligning reads to T7 reference..." >> {log.out}
        bowtie2 -p {threads} --very-sensitive -x results/T7_ref/build_index/T7 \
            -1 {input.r1} -2 {input.r2} \
            -S results/2_T7_removal/1_sam/{wildcards.sample}/{wildcards.sample}_mapped_and_unmapped.sam \
            >> {log.out} 2>> {log.err}
        
        # Convert SAM to BAM
        echo "Step 2: Converting SAM to BAM..." >> {log.out}
        samtools view -@ {threads} -bS results/2_T7_removal/1_sam/{wildcards.sample}/{wildcards.sample}_mapped_and_unmapped.sam > \
            results/2_T7_removal/2_bam/{wildcards.sample}/{wildcards.sample}_mapped_and_unmapped.bam \
            2>> {log.err}
        
        # Extract unmapped reads
        echo "Step 3: Extracting unmapped reads..." >> {log.out}
        samtools view -b -f 12 -F 256 -@ {threads} \
            results/2_T7_removal/2_bam/{wildcards.sample}/{wildcards.sample}_mapped_and_unmapped.bam > \
            results/2_T7_removal/3_unmapped/{wildcards.sample}/{wildcards.sample}_bothReadsUnmapped.bam \
            2>> {log.err}
        
        # Sort unmapped reads
        echo "Step 4: Sorting unmapped reads..." >> {log.out}
        samtools sort -n -m {params.sort_mem} -@ {threads} \
            results/2_T7_removal/3_unmapped/{wildcards.sample}/{wildcards.sample}_bothReadsUnmapped.bam -o \
            results/2_T7_removal/4_unmapped_sorted/{wildcards.sample}/{wildcards.sample}_bothReadsUnmapped_sorted.bam \
            2>> {log.err}
        
        # Convert back to FASTQ
        echo "Step 5: Converting back to FASTQ..." >> {log.out}
        samtools fastq -@ {threads} \
            results/2_T7_removal/4_unmapped_sorted/{wildcards.sample}/{wildcards.sample}_bothReadsUnmapped_sorted.bam \
            -1 {output.r1_clean} -2 {output.r2_clean} \
            -0 /dev/null -s /dev/null -n \
            >> {log.out} 2>> {log.err}
        
        # Clean up intermediate files to save space
        echo "Step 6: Cleaning up intermediate files..." >> {log.out}
        rm -rf results/2_T7_removal/1_sam/{wildcards.sample} \
               results/2_T7_removal/2_bam/{wildcards.sample} \
               results/2_T7_removal/3_unmapped/{wildcards.sample} \
               results/2_T7_removal/4_unmapped_sorted/{wildcards.sample}
        
        echo "T7 removal completed successfully for sample {wildcards.sample}" >> {log.out}
        """
#Step 4: Assembly with SPAdes
rule spades_assembly:
    input:
        r1="results/2_T7_removal/5_removed_sequence/{sample}/{sample}_host_removed_R1.fastq.gz",
        r2="results/2_T7_removal/5_removed_sequence/{sample}/{sample}_host_removed_R2.fastq.gz"
    output:
        scaffolds="results/3_spades_result/{sample}/{sample}_no_correction/scaffolds.fasta"
    conda: "envs/spades.yaml"
    log:
        out="log/4_spades_assembly/{sample}_spades_assembly.log",
        err="log/4_spades_assembly/{sample}_spades_assembly.err"
    threads: config["spades"]["threads"]
    resources:
        slurm_partition=config["high_partition"],
        runtime=config["runtime"],
        mem_mb_per_cpu=config["high_memory"],
        cpus_per_task=config["spades"]["threads"],
        slurm_account=config["slurm_account"]
    params:
        memory=config["spades"]["memory"],
        k_list=config["spades"]["k_values"]
    shell:
        """
        mkdir -p results/3_spades_result/{wildcards.sample}
        mkdir -p log/4_spades_assembly
        spades.py --meta \
            -o results/3_spades_result/{wildcards.sample}/{wildcards.sample}_no_correction \
            -1 {input.r1} -2 {input.r2} \
            -t {threads} -m {params.memory} \
            -k {params.k_list} \
            --only-assembler \
            > {log.out} 2> {log.err}
        """
# Step 5: Rename and filter assemblies
rule rename_filter_assemblies:
    input:
        scaffolds="results/3_spades_result/{sample}/{sample}_no_correction/scaffolds.fasta"
    output:
        renamed="results/4_rename_assembly/rename_5000/{sample}_scaffolds_rename_5000.fasta"
    conda: "envs/seqkit.yaml"
    log:
        out="log/5_rename_filter/{sample}_rename_filter.log",
        err="log/5_rename_filter/{sample}_rename_filter.err"
    threads: config["seqkit"]["threads"]
    resources:
        slurm_partition=config["regular_partition"],
        runtime=config["runtime"],
        mem_mb_per_cpu=config["regular_memory"],
        cpus_per_task=config["seqkit"]["threads"],
        slurm_account=config["slurm_account"]
    params:
        min_length=config["seqkit"]["min_length"],
        nr_width=config["seqkit"]["nr_width"]
    shell:
        """
        mkdir -p results/4_rename_assembly/{{spade_result,rename_5000}}
        mkdir -p log/5_rename_filter


        # Filter sequences >= min_length and rename
        seqkit seq -m {params.min_length} {input.scaffolds} | \
        seqkit replace -p .+ -r "{wildcards.sample}_{{nr}}" --nr-width {params.nr_width} \
        -o {output.renamed} \
        > {log.out} 2> {log.err}
        """
rule checkv_analysis:
    input:
        cluster="results/4_rename_assembly/rename_5000/{sample}_scaffolds_rename_5000.fasta"
    output:
        directory("results/9_checkv/checkv_sample_figure/{sample}")
    conda: "envs/checkv.yaml"
    log:
        out="log/6_checkv/{sample}_checkv.log",
        err="log/6_checkv/{sample}_checkv.err"
    threads: config["checkv"]["threads"]
    resources:
        slurm_partition=config["regular_partition"],
        runtime=config["runtime"],
        mem_mb_per_cpu=config["regular_memory"],
        cpus_per_task=config["checkv"]["threads"],
        slurm_account=config["slurm_account"]
    params:
        db=config["checkv_db"]
    shell:
        """
        mkdir -p results/9_checkv/checkv_sample_figure/{wildcards.sample}
        mkdir -p log/6_checkv
        checkv end_to_end {input.cluster} {output} \
            -d {params.db} \
            -t {threads} \
            > {log.out} 2> {log.err}
        """


rule genomad_identification:
    input:
        cluster="results/4_rename_assembly/rename_5000/{sample}_scaffolds_rename_5000.fasta"
    output:
        directory("results/5_genomad/{sample}")
    conda: "envs/genomad.yaml"
    log:
        out="log/7_genomad/{sample}_genomad.log",
        err="log/7_genomad/{sample}_genomad.err"
    threads: config["genomad"]["threads"]
    resources:
        slurm_partition=config["regular_partition"],
        runtime=config["runtime"],
        mem_mb_per_cpu=config["regular_memory"],
        cpus_per_task=config["genomad"]["threads"],
        slurm_account=config["slurm_account"]
    params:
        db=config["genomad_db"]
    shell:
        """
        mkdir -p results/5_genomad/{wildcards.sample}
        mkdir -p log/7_genomad
        genomad end-to-end --cleanup -t {threads} \
            {input.cluster} \
            {output} \
            {params.db} \
            > {log.out} 2> {log.err}
        """


rule virsorter2_identification:
    input:
        cluster="results/4_rename_assembly/rename_5000/{sample}_scaffolds_rename_5000.fasta"
    output:
        directory("results/6_virsorter2/{sample}")
    conda: "envs/vs2.yaml"
    log:
        out="log/8_virsorter2/{sample}_vs2.log",
        err="log/8_virsorter2/{sample}_vs2.err"
    threads: config["virsorter2"]["threads"]
    resources:
        slurm_partition=config["regular_partition"],
        runtime=config["runtime"],
        mem_mb_per_cpu=config["regular_memory"],
        cpus_per_task=config["virsorter2"]["threads"],
        slurm_account=config["slurm_account"]
    params:
        min_length=config["virsorter2"]["min_length"],
        min_score=config["virsorter2"]["min_score"],
        groups=config["virsorter2"]["groups"],
        database=config["virsorter2"]["database"]
    shell:
        """
        mkdir -p results/6_virsorter2/{wildcards.sample}
        mkdir -p results/6_virsorter2/{wildcards.sample}/tmp
        mkdir -p log/8_virsorter2
        virsorter run --keep-original-seq \
            -i {input.cluster} \
            -w {output} \
            --tmpdir results/6_virsorter2/{wildcards.sample}/tmp \
            -d {params.database} \
            --include-groups {params.groups} \
            --min-length {params.min_length} \
            --min-score {params.min_score} \
            -j {threads} all \
            > {log.out} 2> {log.err}
        """

# Step 9: Viral identification with DeepVirFinder
rule deepvirfinder_identification:
    input:
        cluster="results/4_rename_assembly/rename_5000/{sample}_scaffolds_rename_5000.fasta"
    output:
        result="results/7_deepvirfinder/{sample}/{sample}_scaffolds_rename_5000.fasta_gt3000bp_dvfpred.txt"
    conda: "envs/dvf.yaml"
    log:
        out="log/9_deepvirfinder/{sample}_dvf.log",
        err="log/9_deepvirfinder/{sample}_dvf.err"
    threads: config["deepvirfinder"]["threads"]
    resources:
        slurm_partition=config["regular_partition"],
        runtime=config["runtime"],
        mem_mb_per_cpu=config["regular_memory"],
        cpus_per_task=config["deepvirfinder"]["threads"],
        slurm_account=config["slurm_account"]
    params:
        dvf_path=config["deepvirfinder"]["dvf_path"],
        dvf_models=config["deepvirfinder"]["dvf_models"],
        min_length=config["deepvirfinder"]["min_length"]
    shell:
        """
        mkdir -p results/7_deepvirfinder/{wildcards.sample}
        mkdir -p log/9_deepvirfinder
        python {params.dvf_path} \
            -i {input.cluster} \
            -o results/7_deepvirfinder/{wildcards.sample}/ \
            -l {params.min_length} \
            -c {threads} \
            -m {params.dvf_models} \
            > {log.out} 2> {log.err}
        """







# "--------------------------------------------------------"



# Step 10: Merge viral identification results
rule merge_viral_results:
    input:
        dvf="results/7_deepvirfinder/{sample}/{sample}_scaffolds_rename_5000.fasta_gt3000bp_dvfpred.txt",
        genomad="results/5_genomad/{sample}",
        vs2="results/6_virsorter2/{sample}"
    output:
        csv="results/8_python_merge_filter/csv/{sample}_merged_results.csv",
        list="results/8_python_merge_filter/list/{sample}_merge3_list.txt"
    conda: "envs/python.yaml"
    log:
        out="log/10_merge_viral/{sample}_merge.log",
        err="log/10_merge_viral/{sample}_merge.err"
    threads: config["merge_filter"]["threads"]
    resources:
        slurm_partition=config["regular_partition"],
        runtime=config["runtime"],
        mem_mb_per_cpu=config["regular_memory"],
        cpus_per_task=config["merge_filter"]["threads"],
        slurm_account=config["slurm_account"]
    params:
        merge_script=config["merge_script"]
    shell:
        """
        mkdir -p results/8_python_merge_filter
        mkdir -p log/10_merge_viral

        # Copy GeNomad summary file to expected location
        mkdir -p results/5_genomad/genomad_result_summary
        cp {input.genomad}/{wildcards.sample}_scaffolds_rename_5000_summary/{wildcards.sample}_scaffolds_rename_5000_virus_summary.tsv \
           results/5_genomad/genomad_result_summary/

        python {params.merge_script} \
            --dvf {input.dvf} \
            --genomad results/5_genomad/genomad_result_summary/{wildcards.sample}_scaffolds_rename_5000_virus_summary.tsv \
            --vs2 {input.vs2}/final-viral-score.tsv \
            --output results/8_python_merge_filter \
            > {log.out} 2> {log.err}
        """


# Step 11: Extract viral sequences
rule extract_viral_sequences:
    input:
        list="results/8_python_merge_filter/list/{sample}_merge3_list.txt",
        fasta="results/4_rename_assembly/rename_5000/{sample}_scaffolds_rename_5000.fasta"
    output:
        filtered="results/8_python_merge_filter/seqkit_filter_vircontig/{sample}_filterd_vircontig.fasta"
    conda: "envs/seqkit.yaml"
    log:
        out="log/11_extract_viral/{sample}_extract.log",
        err="log/11_extract_viral/{sample}_extract.err"
    threads: config["seqkit"]["threads"]
    resources:
        slurm_partition=config["regular_partition"],
        runtime=config["runtime"],
        mem_mb_per_cpu=config["regular_memory"],
        cpus_per_task=config["seqkit"]["threads"],
        slurm_account=config["slurm_account"]
    shell:
        """
        mkdir -p results/8_python_merge_filter/seqkit_filter_vircontig
        mkdir -p log/11_extract_viral
        seqkit grep -f {input.list} {input.fasta} > {output.filtered} \
            2> {log.err}
        """

# Step 12: Second CheckV analysis on filtered sequences
rule checkv_filtered_analysis:
    input:
        filtered="results/8_python_merge_filter/seqkit_filter_vircontig/{sample}_filterd_vircontig.fasta"
    output:
        summary="results/9_checkv/chekv_result/{sample}/quality_summary.tsv"
    conda: "envs/checkv.yaml"
    log:
        out="log/12_checkv_filtered/{sample}_checkv_filtered.log",
        err="log/12_checkv_filtered/{sample}_checkv_filtered.err"
    threads: config["checkv"]["threads"]
    resources:
        slurm_partition=config["regular_partition"],
        runtime=config["runtime"],
        mem_mb_per_cpu=config["regular_memory"],
        cpus_per_task=config["checkv"]["threads"],
        slurm_account=config["slurm_account"]
    params:
        db=config["checkv_db"]
    shell:
        """
        mkdir -p results/9_checkv/chekv_result/{wildcards.sample}
        mkdir -p log/12_checkv_filtered
        checkv end_to_end {input.filtered} results/9_checkv/chekv_result/{wildcards.sample} \
            -d {params.db} \
            -t {threads} \
            > {log.out} 2> {log.err}
        """

# Step 13: Python filter based on CheckV results
rule python_checkv_filter:
    input:
        genomad="results/5_genomad/{sample}",
        vs2="results/6_virsorter2/{sample}",
        checkv="results/9_checkv/chekv_result/{sample}/quality_summary.tsv"
    output:
        list="results/9_checkv/filtered_list/{sample}_checkv_extract.txt"
    conda: "envs/python.yaml"
    log:
        out="log/13_python_filter/{sample}_python_filter.log",
        err="log/13_python_filter/{sample}_python_filter.err"
    threads: config["merge_filter"]["threads"]
    resources:
        slurm_partition=config["regular_partition"],
        runtime=config["runtime"],
        mem_mb_per_cpu=config["regular_memory"],
        cpus_per_task=config["merge_filter"]["threads"],
        slurm_account=config["slurm_account"]
    params:
        checkv_filter_script=config["checkv_filter_script"]
    shell:
        """
        mkdir -p results/9_checkv/filtered_list
        mkdir -p log/13_python_filter
        python {params.checkv_filter_script} \
            --checkv {input.checkv} \
            --genomad results/5_genomad/genomad_result_summary/{wildcards.sample}_scaffolds_rename_5000_virus_summary.tsv \
            --vs2 {input.vs2}/final-viral-score.tsv \
            --output results/9_checkv/filtered_list \
            > {log.out} 2> {log.err}
        touch {output.list}
        """

# Step 14: Extract final filtered viral sequences
rule extract_final_viral_sequences:
    input:
        list="results/9_checkv/filtered_list/{sample}_checkv_extract.txt",
        fasta="results/4_rename_assembly/rename_5000/{sample}_scaffolds_rename_5000.fasta"
    output:
        final="results/10_cluster/1_filtered_seq/{sample}_filterd_final_vircontig.fasta"
    conda: "envs/seqkit.yaml"
    log:
        out="log/14_extract_final/{sample}_extract_final.log",
        err="log/14_extract_final/{sample}_extract_final.err"
    threads: config["seqkit"]["threads"]
    resources:
        slurm_partition=config["regular_partition"],
        runtime=config["runtime"],
        mem_mb_per_cpu=config["regular_memory"],
        cpus_per_task=config["seqkit"]["threads"],
        slurm_account=config["slurm_account"]
    shell:
        """
        mkdir -p results/10_cluster/1_filtered_seq
        mkdir -p log/14_extract_final
        seqkit grep -f {input.list} {input.fasta} > {output.final} \
            2> {log.err}
        """


# Step 15a: Cluster individual samples using CheckV method
rule cluster_per_sample:
    input:
        fasta="results/10_cluster/1_filtered_seq/{sample}_filterd_final_vircontig.fasta"
    output:
        blast_db=directory("results/10_cluster/3_clustered/per_sample/{sample}/blast_db"),
        blast_results="results/10_cluster/3_clustered/per_sample/{sample}/{sample}_blast.tsv",
        ani_results="results/10_cluster/3_clustered/per_sample/{sample}/{sample}_ani.tsv",
        cluster_results="results/10_cluster/3_clustered/per_sample/{sample}/{sample}_cluster.tsv"
    conda: "envs/checkv.yaml"
    threads: 16
    log:
        out="log/15_cluster_per_sample/{sample}_cluster.log",
        err="log/15_cluster_per_sample/{sample}_cluster.err"
    resources:
        slurm_partition="compute,memory",
        runtime = 7200,
        mem_mb_per_cpu = 3900,
        cpus_per_task = 32,
        slurm_account = "research-as-bt"
    shell:
        """
        mkdir -p results/10_cluster/3_clustered/per_sample/{wildcards.sample}
        mkdir -p log/15_cluster_per_sample

        # Create BLAST database
        makeblastdb -in {input.fasta} \
            -dbtype nucl \
            -out results/10_cluster/3_clustered/per_sample/{wildcards.sample}/{wildcards.sample}_db \
            >> {log.out} 2>> {log.err}

        # Run BLAST all-vs-all
        blastn -query {input.fasta} \
            -db results/10_cluster/3_clustered/per_sample/{wildcards.sample}/{wildcards.sample}_db \
            -outfmt '6 std qlen slen' \
            -max_target_seqs 10000 \
            -out {output.blast_results} \
            -num_threads {threads} \
            >> {log.out} 2>> {log.err}

        # Calculate ANI
        python {config[checkv_scripts]}/anicalc.py \
            -i {output.blast_results} \
            -o {output.ani_results} \
            >> {log.out} 2>> {log.err}

        # Cluster sequences
        python {config[checkv_scripts]}/aniclust.py \
            --fna {input.fasta} \
            --ani {output.ani_results} \
            --out {output.cluster_results} \
            --min_ani 95 \
            --min_tcov 85 \
            --min_qcov 0 \
            >> {log.out} 2>> {log.err}

        # Create directory for blast database files and move them
        mkdir -p {output.blast_db}
        mv results/10_cluster/3_clustered/per_sample/{wildcards.sample}/{wildcards.sample}_db.* {output.blast_db}/
        """

# Step 15b: Cluster all samples using CheckV method
rule cluster_all_samples:
    input:
        expand("results/10_cluster/1_filtered_seq/{sample}_filterd_final_vircontig.fasta", sample=config["samples"])
    output:
        combined="results/10_cluster/2_combined/all_samples_combined.fasta",
        blast_db=directory("results/10_cluster/3_clustered/blast_db"),
        blast_results="results/10_cluster/3_clustered/all_samples_blast.tsv",
        ani_results="results/10_cluster/3_clustered/all_samples_ani.tsv",
        cluster_results="results/10_cluster/3_clustered/all_samples_cluster.tsv"
    conda: "envs/checkv.yaml"
    threads: 16
    log:
        out="log/15_cluster_all/cluster_all.log",
        err="log/15_cluster_all/cluster_all.err"
    resources:
        slurm_partition="compute,memory",
        runtime = 7200,
        mem_mb_per_cpu = 3900,
        cpus_per_task = 32,
        slurm_account = "research-as-bt"
    shell:
        """
        mkdir -p results/10_cluster/2_combined
        mkdir -p results/10_cluster/3_clustered
        mkdir -p log/15_cluster_all

        # Combine all filtered sequences
        cat {input} > {output.combined}

        # Create BLAST database
        makeblastdb -in {output.combined} \
            -dbtype nucl \
            -out results/10_cluster/3_clustered/all_samples_db \
            >> {log.out} 2>> {log.err}

        # Run BLAST all-vs-all
        blastn -query {output.combined} \
            -db results/10_cluster/3_clustered/all_samples_db \
            -outfmt '6 std qlen slen' \
            -max_target_seqs 10000 \
            -out {output.blast_results} \
            -num_threads {threads} \
            >> {log.out} 2>> {log.err}

        # Calculate ANI
        python {config[checkv_scripts]}/anicalc.py \
            -i {output.blast_results} \
            -o {output.ani_results} \
            >> {log.out} 2>> {log.err}

        # Cluster sequences
        python {config[checkv_scripts]}/aniclust.py \
            --fna {output.combined} \
            --ani {output.ani_results} \
            --out {output.cluster_results} \
            --min_ani 95 \
            --min_tcov 85 \
            --min_qcov 0 \
            >> {log.out} 2>> {log.err}

        # Create directory for blast database files
        mkdir -p {output.blast_db}
        mv results/10_cluster/3_clustered/all_samples_db.* {output.blast_db}/
        """

# Step 16a: Extract representative sequences per sample
rule extract_representatives_per_sample:
    input:
        fasta="results/10_cluster/1_filtered_seq/{sample}_filterd_final_vircontig.fasta",
        cluster_results="results/10_cluster/3_clustered/per_sample/{sample}/{sample}_cluster.tsv"
    output:
        representatives="results/10_cluster/3_clustered/per_sample/{sample}/{sample}_cluster_representatives.fasta",
        copied="results/10_cluster/4_cluster_fasta/{sample}_cluster_representatives.fasta"
    conda: "envs/seqkit.yaml"
    log:
        out="log/16_extract_representatives/{sample}_extract_representatives.log",
        err="log/16_extract_representatives/{sample}_extract_representatives.err"
    resources:
        slurm_partition="compute",
        runtime=240,
        mem_mb_per_cpu=3000,
        cpus_per_task=2,
        slurm_account="research-as-bt"
    shell:
        """
        mkdir -p log/16_extract_representatives
        mkdir -p results/10_cluster/4_cluster_fasta

        # Extract representative sequences (cluster centroids)
        cut -f 2 {input.cluster_results} | seqkit grep -f - -w 0 {input.fasta} > {output.representatives} \
            2> {log.err}

        # Copy to 4_cluster_fasta directory
        cp {output.representatives} {output.copied}
        """

# Step 16b: Extract representative sequences for all samples
rule extract_representatives_all_samples:
    input:
        combined="results/10_cluster/2_combined/all_samples_combined.fasta",
        cluster_results="results/10_cluster/3_clustered/all_samples_cluster.tsv"
    output:
        representatives="results/10_cluster/3_clustered/all_samples_cluster_representatives.fasta",
        copied="results/10_cluster/4_cluster_fasta/all_samples_cluster_representatives.fasta"
    conda: "envs/seqkit.yaml"
    log:
        out="log/16_extract_representatives/all_samples_extract_representatives.log",
        err="log/16_extract_representatives/all_samples_extract_representatives.err"
    resources:
        slurm_partition="compute",
        runtime=240,
        mem_mb_per_cpu=3000,
        cpus_per_task=2,
        slurm_account="research-as-bt"
    shell:
        """
        mkdir -p log/16_extract_representatives
        mkdir -p results/10_cluster/4_cluster_fasta

        # Extract representative sequences (cluster centroids)
        cut -f 2 {input.cluster_results} | seqkit grep -f - -w 0 {input.combined} > {output.representatives} \
            2> {log.err}

        # Copy to 4_cluster_fasta directory
        cp {output.representatives} {output.copied}
        """

# Step 17a: Build bowtie2 index for per-sample representative sequences
rule build_cluster_index_per_sample:
    input:
        representatives="results/10_cluster/4_cluster_fasta/{sample}_cluster_representatives.fasta"
    output:
        flag="results/10_cluster/5_bowtie_index/per_sample/{sample}/cluster_index.done"
    conda: "envs/bowtie2.yaml"
    log:
        out="log/17_build_cluster_index/{sample}_build_index.log",
        err="log/17_build_cluster_index/{sample}_build_index.err"
    threads: config["cluster_index_build"]["threads"]
    resources:
        slurm_partition=config["cluster_index_build"]["partition"],
        runtime=config["cluster_index_build"]["runtime"],
        mem_mb_per_cpu=config["cluster_index_build"]["mem_mb_per_cpu"],
        cpus_per_task=config["cluster_index_build"]["threads"],
        slurm_account=config["slurm_account"]
    shell:
        """
        mkdir -p results/10_cluster/5_bowtie_index/per_sample/{wildcards.sample}
        mkdir -p log/17_build_cluster_index

        # Build bowtie2 index
        bowtie2-build {input.representatives} results/10_cluster/5_bowtie_index/per_sample/{wildcards.sample}/cluster_index \
            > {log.out} 2> {log.err}

        touch {output.flag}
        """

# Step 17b: Build bowtie2 index for all samples representative sequences
rule build_cluster_index_all_samples:
    input:
        representatives="results/10_cluster/4_cluster_fasta/all_samples_cluster_representatives.fasta"
    output:
        flag="results/10_cluster/5_bowtie_index/all_samples/cluster_index.done"
    conda: "envs/bowtie2.yaml"
    log:
        out="log/17_build_cluster_index/all_samples_build_index.log",
        err="log/17_build_cluster_index/all_samples_build_index.err"
    threads: config["cluster_index_build"]["threads"]
    resources:
        slurm_partition=config["cluster_index_build"]["partition"],
        runtime=config["cluster_index_build"]["runtime"],
        mem_mb_per_cpu=config["cluster_index_build"]["mem_mb_per_cpu"],
        cpus_per_task=config["cluster_index_build"]["threads"],
        slurm_account=config["slurm_account"]
    shell:
        """
        mkdir -p results/10_cluster/5_bowtie_index/all_samples
        mkdir -p log/17_build_cluster_index

        # Build bowtie2 index
        bowtie2-build {input.representatives} results/10_cluster/5_bowtie_index/all_samples/cluster_index \
            > {log.out} 2> {log.err}

        touch {output.flag}
        """

# Step 18a: Bowtie2 alignment for per-sample clusters
rule bowtie2_alignment_per_sample:
    input:
        r1="results/2_T7_removal/5_removed_sequence/{sample}/{sample}_host_removed_R1.fastq.gz",
        r2="results/2_T7_removal/5_removed_sequence/{sample}/{sample}_host_removed_R2.fastq.gz",
        index_flag="results/10_cluster/5_bowtie_index/per_sample/{sample}/cluster_index.done"
    output:
        bam="results/11_bowtie2/per_sample/{sample}/{sample}.f2.sorted.bam",
        bam_index="results/11_bowtie2/per_sample/{sample}/{sample}.f2.sorted.bam.bai"
    conda: "envs/bowtie2_samtools.yaml"
    threads: config["bowtie2_alignment"]["threads"]
    log:
        out="log/18_bowtie2_alignment/per_sample/{sample}_alignment.log",
        err="log/18_bowtie2_alignment/per_sample/{sample}_alignment.err"
    resources:
        slurm_partition=config["bowtie2_alignment"]["partition"],
        runtime=config["bowtie2_alignment"]["runtime"],
        mem_mb_per_cpu=config["bowtie2_alignment"]["mem_mb_per_cpu"],
        cpus_per_task=config["bowtie2_alignment"]["threads"],
        slurm_account=config["slurm_account"]
    shell:
        """
        mkdir -p results/11_bowtie2/per_sample/{wildcards.sample}
        mkdir -p log/18_bowtie2_alignment/per_sample

        # Bowtie2 alignment
        bowtie2 --sensitive -t -p {threads} \
            -x results/10_cluster/5_bowtie_index/per_sample/{wildcards.sample}/cluster_index \
            -1 {input.r1} -2 {input.r2} \
            -S results/11_bowtie2/per_sample/{wildcards.sample}/{wildcards.sample}.sam \
            > {log.out} 2> {log.err}

        # Samtools processing
        samtools view -@ {threads} -hbS -f 2 results/11_bowtie2/per_sample/{wildcards.sample}/{wildcards.sample}.sam \
        | samtools sort -@ {threads} -o {output.bam} - \
            >> {log.out} 2>> {log.err}

        samtools index -@ {threads} {output.bam} \
            >> {log.out} 2>> {log.err}

        # Clean up SAM file to save space
        rm -f results/11_bowtie2/per_sample/{wildcards.sample}/{wildcards.sample}.sam
        """

# Step 18b: Bowtie2 alignment for all samples clusters (all samples align to combined clusters)
rule bowtie2_alignment_all_samples:
    input:
        r1="results/2_T7_removal/5_removed_sequence/{sample}/{sample}_host_removed_R1.fastq.gz",
        r2="results/2_T7_removal/5_removed_sequence/{sample}/{sample}_host_removed_R2.fastq.gz",
        index_flag="results/10_cluster/5_bowtie_index/all_samples/cluster_index.done"
    output:
        bam="results/11_bowtie2/all_samples/{sample}/{sample}.f2.sorted.bam",
        bam_index="results/11_bowtie2/all_samples/{sample}/{sample}.f2.sorted.bam.bai"
    conda: "envs/bowtie2_samtools.yaml"
    threads: config["bowtie2_alignment"]["threads"]
    log:
        out="log/18_bowtie2_alignment/all_samples/{sample}_alignment.log",
        err="log/18_bowtie2_alignment/all_samples/{sample}_alignment.err"
    resources:
        slurm_partition=config["bowtie2_alignment"]["partition"],
        runtime=config["bowtie2_alignment"]["runtime"],
        mem_mb_per_cpu=config["bowtie2_alignment"]["mem_mb_per_cpu"],
        cpus_per_task=config["bowtie2_alignment"]["threads"],
        slurm_account=config["slurm_account"]
    shell:
        """
        mkdir -p results/11_bowtie2/all_samples/{wildcards.sample}
        mkdir -p log/18_bowtie2_alignment/all_samples

        # Bowtie2 alignment
        bowtie2 --sensitive -t -p {threads} \
            -x results/10_cluster/5_bowtie_index/all_samples/cluster_index \
            -1 {input.r1} -2 {input.r2} \
            -S results/11_bowtie2/all_samples/{wildcards.sample}/{wildcards.sample}.sam \
            > {log.out} 2> {log.err}

        # Samtools processing
        samtools view -@ {threads} -hbS -f 2 results/11_bowtie2/all_samples/{wildcards.sample}/{wildcards.sample}.sam \
        | samtools sort -@ {threads} -o {output.bam} - \
            >> {log.out} 2>> {log.err}

        samtools index -@ {threads} {output.bam} \
            >> {log.out} 2>> {log.err}

        # Clean up SAM file to save space
        rm -f results/11_bowtie2/all_samples/{wildcards.sample}/{wildcards.sample}.sam
        """

# Step 19a: Calculate abundance for per-sample clusters
rule calculate_abundance_per_sample:
    input:
        bam="results/11_bowtie2/per_sample/{sample}/{sample}.f2.sorted.bam",
        bam_index="results/11_bowtie2/per_sample/{sample}/{sample}.f2.sorted.bam.bai"
    output:
        tpm="results/11_bowtie2/per_sample/{sample}/{sample}_TPM.tsv",
        mean="results/11_bowtie2/per_sample/{sample}/{sample}_mean.tsv",
        read_count="results/11_bowtie2/per_sample/{sample}/{sample}_count.tsv"
    conda: "envs/coverm.yaml"
    threads: config["coverm"]["threads"]
    log:
        out="log/19_coverm/per_sample/{sample}_coverm.log",
        err="log/19_coverm/per_sample/{sample}_coverm.err"
    resources:
        slurm_partition=config["coverm"]["partition"],
        runtime=config["coverm"]["runtime"],
        mem_mb_per_cpu=config["coverm"]["mem_mb_per_cpu"],
        cpus_per_task=config["coverm"]["threads"],
        slurm_account=config["slurm_account"]
    shell:
        """
        mkdir -p log/19_coverm/per_sample

        # CoverM abundance calculations
        coverm contig -m tpm \
            --bam-files {input.bam} \
            --min-read-percent-identity 95 \
            --min-read-aligned-percent 90 \
            --output-file {output.tpm} \
            --contig-end-exclusion 0 \
            --no-zeros -t {threads} \
            > {log.out} 2> {log.err}

        coverm contig -m mean \
            --bam-files {input.bam} \
            --min-read-percent-identity 95 \
            --min-read-aligned-percent 90 \
            --output-file {output.mean} \
            --contig-end-exclusion 0 \
            --no-zeros -t {threads} \
            >> {log.out} 2>> {log.err}

        coverm contig -m count \
            --bam-files {input.bam} \
            --min-read-percent-identity 95 \
            --min-read-aligned-percent 90 \
            --output-file {output.read_count} \
            --contig-end-exclusion 0 \
            --no-zeros -t {threads} \
            >> {log.out} 2>> {log.err}
        """

# Step 19b: Calculate abundance for all samples clusters
rule calculate_abundance_all_samples:
    input:
        bam="results/11_bowtie2/all_samples/{sample}/{sample}.f2.sorted.bam",
        bam_index="results/11_bowtie2/all_samples/{sample}/{sample}.f2.sorted.bam.bai"
    output:
        tpm="results/11_bowtie2/all_samples/{sample}/{sample}_TPM.tsv",
        mean="results/11_bowtie2/all_samples/{sample}/{sample}_mean.tsv",
        read_count="results/11_bowtie2/all_samples/{sample}/{sample}_count.tsv"
    conda: "envs/coverm.yaml"
    threads: config["coverm"]["threads"]
    log:
        out="log/19_coverm/all_samples/{sample}_coverm.log",
        err="log/19_coverm/all_samples/{sample}_coverm.err"
    resources:
        slurm_partition=config["coverm"]["partition"],
        runtime=config["coverm"]["runtime"],
        mem_mb_per_cpu=config["coverm"]["mem_mb_per_cpu"],
        cpus_per_task=config["coverm"]["threads"],
        slurm_account=config["slurm_account"]
    shell:
        """
        mkdir -p log/19_coverm/all_samples

        # CoverM abundance calculations
        coverm contig -m tpm \
            --bam-files {input.bam} \
            --min-read-percent-identity 95 \
            --min-read-aligned-percent 90 \
            --output-file {output.tpm} \
            --contig-end-exclusion 0 \
            --no-zeros -t {threads} \
            > {log.out} 2> {log.err}

        coverm contig -m mean \
            --bam-files {input.bam} \
            --min-read-percent-identity 95 \
            --min-read-aligned-percent 90 \
            --output-file {output.mean} \
            --contig-end-exclusion 0 \
            --no-zeros -t {threads} \
            >> {log.out} 2>> {log.err}

        coverm contig -m count \
            --bam-files {input.bam} \
            --min-read-percent-identity 95 \
            --min-read-aligned-percent 90 \
            --output-file {output.read_count} \
            --contig-end-exclusion 0 \
            --no-zeros -t {threads} \
            >> {log.out} 2>> {log.err}
        """
