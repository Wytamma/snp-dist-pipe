import os

# Reference genome
REF_GENOME = config.get("reference")

# Mask file (update the path to your mask file)
MASK_FILE = config.get('mask')

EXTENSION = config.get('extension')

# Path to the existing alignment of ~44k samples
EXISTING_ALIGNMENT = config.get('alignment')
if not EXISTING_ALIGNMENT:
    raise ValueError("An existing alignment is required!")
# Assemblies
ASSEMBLY_DIR = config.get('assemblies')
if not ASSEMBLY_DIR:
    raise ValueError("Assemblies are required!")

# Get sample names by parsing the read filenames
SAMPLES, = glob_wildcards(os.path.join(ASSEMBLY_DIR, "{sample}" + EXTENSION))

print(f"Found {len(SAMPLES)} samples in {ASSEMBLY_DIR} with extension {EXTENSION}")
print(f"REF_GENOME: {REF_GENOME}")
print(f"MASK_FILE: {MASK_FILE}")
print(f"EXISTING_ALIGNMENT: {EXISTING_ALIGNMENT}")
print(f"ASSEMBLY_DIR: {ASSEMBLY_DIR}")
print(f"EXTENSION: {EXTENSION}")

rule all:
    input:
        tidy_distance_matrix = "tidy_distance_matrix.tsv"

rule snippy:
    input:
        assembly = os.path.join(ASSEMBLY_DIR, "{sample}" + EXTENSION),
        ref = REF_GENOME
    output:
        vcf = temp("snippy/{sample}/snps.vcf.gz")  # temp() marks it as intermediate
    params:
        outdir = "snippy/{sample}"
    conda:
        "envs/snippy.yaml"
    threads: 4
    shadow: "minimal"
    shell:
        """
        snippy --outdir {params.outdir} --ref {input.ref} --ctgs {input.assembly} --cpus {threads} --force
        """

rule bcftools_consensus:
    input:
        ref = REF_GENOME,
        vcf = "snippy/{sample}/snps.vcf.gz"
    output:
        consensus = temp("consensus/{sample}.fa"),  # temp() marks it as intermediate
        filtered_snps = temp("consensus/{sample}.filtered_snps.vcf.gz"),
        filtered_snps_index = temp("consensus/{sample}.filtered_snps.vcf.gz.csi")
    conda:
        "envs/bcftools.yaml"
    shadow: "minimal"
    shell:
        """
        # First, filter the VCF to include only SNPs
        bcftools view -v snps {input.vcf} -Oz -o {output.filtered_snps}
        
        # Index the filtered VCF file
        bcftools index {output.filtered_snps}

        # Now, apply the filtered and indexed VCF to the reference genome
        bcftools consensus -f {input.ref} {output.filtered_snps} > {output.consensus}
        """

# Rule to mask consensus sequences using BEDTools
rule mask_consensus:
    input:
        consensus = "consensus/{sample}.fa",
        mask = MASK_FILE
    output:
        temp_masked = temp("masked/{sample}.temp"),  # temp() marks it as intermediate
        masked = temp("masked/{sample}.fa")  # temp() marks it as intermediate
    conda:
        "envs/bedtools.yaml"
    shadow: "minimal"
    shell:
        """
        bedtools maskfasta -fi {input.consensus} -bed {input.mask} -fo {output.temp_masked}
        sed 's/^>.*/>{wildcards.sample}/' {output.temp_masked} > {output.masked}
        """

# Rule to concatenate all masked consensus sequences into a single MSA
rule concatenate_alignments:
    input:
        masked_sequences = expand("masked/{sample}.fa", sample=SAMPLES)
    output:
        msa = temp("alignment.fa")  # temp() marks it as intermediate
    shadow: "minimal"
    run:
        with open(output.msa, 'w') as outfile:
            for infile in input.masked_sequences:
                with open(infile) as seqfile:
                    outfile.write(seqfile.read())

# Rule to calculate SNP distances using psdm
rule calculate_snp_distances:
    input:
        msa = "alignment.fa",
        existing_msa = EXISTING_ALIGNMENT
    output:
        distance_matrix = temp("distance_matrix.tsv")  # temp() marks it as intermediate
    conda:
        "envs/psdm.yaml"
    threads: workflow.cores
    shadow: "minimal"
    shell:
        """
        psdm -P -t {threads} {input.msa} {input.existing_msa} > {output.distance_matrix}
        """

rule tidy_distance_matrix:
    input:
        distance_matrix = "distance_matrix.tsv"
    output:
        tidy_distance_matrix = "tidy_distance_matrix.tsv"
    conda:
        "envs/python.yaml"
    script: """scripts/tidy_distance_matrix.py"""
