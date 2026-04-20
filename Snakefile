rule all:
    input:
        auspice_json = "auspice/MTB-tree-2.json"
        #auspice_json = "auspice/MTB.json",

input_genotype = "data/genotype_new.tsv",
colors = "config/color.tsv",
auspice_config = "config/auspice_config.json"

rule pull_sequences:
    message:
        """
        Pulling fasta sequences...
        """
    input:
        genotype = input_genotype,
    output:
        multifasta = "data/sequences.fasta",
        ref = "data/reference.fasta"
    shell:
        """
        python src/pull_fasta.py {input.genotype} {output.multifasta} {output.ref}
        """
rule pull_annotations:
    message:
        """
        Pulling sequence annotations...
        """
    input:
        genotype = input_genotype,
    output:
        gbk = "data/sequences.gbk",
        ref_gbk = "data/reference.gbk"
    shell:
        """
        python src/pull_annotations.py {input.genotype} {output.gbk} {output.ref_gbk}
        """
rule extract_metadata:
    message:
        """
        Extracting metadata
        """
    input:
        gbk = rules.pull_annotations.output.gbk,
        genotype = input_genotype
    output:
        metadata = "data/metadata.tsv",

    shell:
        """
        python src/extract_metadata.py {input.gbk} {input.genotype} {output.metadata}
        """

rule align:
    message:
        """
        Aligning sequences to {input.reference}
        """
    input:
        sequences = rules.pull_sequences.output.multifasta,
        reference = rules.pull_sequences.output.ref
    output:
        alignment = "results/aligned.sam"
    conda:
        "minimap"
    shell:
        """
        minimap2 -ax asm5 {input.reference} {input.sequences} > {output}
        """
rule sam2fasta:
    message:
        """
        Converting sam format alignment to multi-sequence alignment format.
        """
    input:
        sam = rules.align.output.alignment,
    output:
        msa = "results/aligned.fasta"
    conda:
        "gofasta"
    shell:
        """
        gofasta sam toMultiAlign -s {input.sam} -o {output.msa}
        """
rule tree:
    message: "Building tree"
    input:
        msa = rules.sam2fasta.output.msa
    output:
        tree= "results/raw_tree.nwk"
    conda: 
        "augur"
    shell:
        """
        augur tree \
        --alignment {input.msa} \
        --output {output.tree}
        """
rule refine:
    message:
        """
        Refining tree
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
          - filter tips more than {params.clock_filter_iqd} IQDs from clock expectation
        """
    input:
        tree = rules.tree.output.tree,
        alignment = rules.sam2fasta.output.msa,
        metadata = rules.extract_metadata.output.metadata
    output:
        tree = "results/tree.nwk",
        node_data = "results/branch_lengths.json"
    params:
        coalescent = "opt",
        date_inference = "marginal",
        clock_filter_iqd = 4
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --alignment {input.alignment} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --coalescent {params.coalescent} \
            --date-confidence \
            --date-inference {params.date_inference} \
            --clock-filter-iqd {params.clock_filter_iqd}
        """
rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine.output.tree,
        alignment = rules.sam2fasta.output.msa
    output:
        node_data = "results/nt_muts.json"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference}
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = rules.pull_annotations.output.ref_gbk
    output:
        node_data = "results/aa_muts.json"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output-node-data {output.node_data} \
        """

rule traits:
    message: "Inferring ancestral traits for {params.columns!s}"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.extract_metadata.output.metadata
    output:
        node_data = "results/traits.json",
    params:
        columns = "clade tbprofiler_lineage"
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output-node-data {output.node_data} \
            --columns {params.columns} \
            --confidence
        """


rule export:
    message: "Exporting data files for Auspice with clade information"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.extract_metadata.output.metadata,
        branch_lengths = rules.refine.output.node_data,
        traits = rules.traits.output.node_data,  # Include traits with clade info
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
        auspice_config = auspice_config
    output:
        auspice_json = rules.all.input.auspice_json,
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.traits} {input.nt_muts} {input.aa_muts} \
            --auspice-config {input.auspice_config} \
            --output {output.auspice_json}
        """
