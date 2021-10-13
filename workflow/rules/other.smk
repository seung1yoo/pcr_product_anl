# An example collection of Snakemake rules imported in the main Snakefile.

def get_fq1(wildcards):
    #print(units.loc['Library_PCR_for_NGS-1'].loc['lane1']['fq1'])
    unit = units.loc[wildcards.sample].loc[wildcards.unit]
    return unit['fq1']

def get_fq2(wildcards):
    unit = units.loc[wildcards.sample].loc[wildcards.unit]
    return unit['fq2']


rule run_flash:
    input:
        fq1 = lambda wildcards: units.loc[wildcards.sample].loc[wildcards.unit, 'fq1'],
        fq2 = lambda wildcards: units.loc[wildcards.sample].loc[wildcards.unit, 'fq2']
    output:
        extended = 'analysis/flash/{sample}_{unit}/{sample}_{unit}.extendedFrags.fastq.gz',
        nc1 = 'analysis/flash/{sample}_{unit}/{sample}_{unit}.notCombined_1.fastq.gz',
        nc2 = 'analysis/flash/{sample}_{unit}/{sample}_{unit}.notCombined_2.fastq.gz',
        hist = 'analysis/flash/{sample}_{unit}/{sample}_{unit}.hist',
        histogram = 'analysis/flash/{sample}_{unit}/{sample}_{unit}.histogram',
    params:
        outdir="analysis/flash/{sample}_{unit}/",
        prefix="{sample}_{unit}",
    threads: 10
    shell:
        "flash"
        " -d {params.outdir}"
        " -o {params.prefix}"
        " --compress"
        " --compress-prog=gzip"
        " --suffix=gz"
        " --threads {threads}"
        " {input.fq1}"
        " {input.fq2}"
        " > {params.outdir}/{params.prefix}.mylog"

rule take_product_ext:
    input:
        fqgz = 'analysis/flash/{sample}_{unit}/{sample}_{unit}.extendedFrags.fastq.gz',
    output:
        fna = 'analysis/product/{sample}_{unit}/{sample}_{unit}.ext.product.fna',
        faa = 'analysis/product/{sample}_{unit}/{sample}_{unit}.ext.product.faa',
        tsv = 'analysis/product/{sample}_{unit}/{sample}_{unit}.ext.product.tsv'
    params:
        primer5 = 'CAAGTGGCCACAAACCACCAGAGT',
        primer3 = 'GCGCAGACCGGCTG',
        prefix = 'analysis/product/{sample}_{unit}/{sample}_{unit}.ext'
    shell:
        'python {config[script_home]}/take_product.py'
        ' --primer5 {params.primer5}'
        ' --primer3 {params.primer3}'
        ' --infqgz {input.fqgz}'
        ' --prefix {params.prefix}'
        
rule take_product_nc1:
    input:
        fqgz = 'analysis/flash/{sample}_{unit}/{sample}_{unit}.notCombined_1.fastq.gz',
    output:
        fna = 'analysis/product/{sample}_{unit}/{sample}_{unit}.nc1.product.fna',
        faa = 'analysis/product/{sample}_{unit}/{sample}_{unit}.nc1.product.faa',
        tsv = 'analysis/product/{sample}_{unit}/{sample}_{unit}.nc1.product.tsv'
    params:
        primer5 = 'CAAGTGGCCACAAACCACCAGAGT',
        primer3 = 'GCGCAGACCGGCTG',
        prefix = 'analysis/product/{sample}_{unit}/{sample}_{unit}.nc1'
    shell:
        'python {config[script_home]}/take_product.py'
        ' --primer5 {params.primer5}'
        ' --primer3 {params.primer3}'
        ' --infqgz {input.fqgz}'
        ' --prefix {params.prefix}'
        
rule take_product_nc2:
    input:
        fqgz = 'analysis/flash/{sample}_{unit}/{sample}_{unit}.notCombined_2.fastq.gz',
    output:
        fna = 'analysis/product/{sample}_{unit}/{sample}_{unit}.nc2.product.fna',
        faa = 'analysis/product/{sample}_{unit}/{sample}_{unit}.nc2.product.faa',
        tsv = 'analysis/product/{sample}_{unit}/{sample}_{unit}.nc2.product.tsv'
    params:
        primer5 = 'CAAGTGGCCACAAACCACCAGAGT',
        primer3 = 'GCGCAGACCGGCTG',
        prefix = 'analysis/product/{sample}_{unit}/{sample}_{unit}.nc2'
    shell:
        'python {config[script_home]}/take_product.py'
        ' --primer5 {params.primer5}'
        ' --primer3 {params.primer3}'
        ' --infqgz {input.fqgz}'
        ' --prefix {params.prefix}'
        
rule merge_fna_product:
    input:
        ext_fna = 'analysis/product/{sample}_{unit}/{sample}_{unit}.ext.product.fna',
        nc1_fna = 'analysis/product/{sample}_{unit}/{sample}_{unit}.nc1.product.fna',
        nc2_fna = 'analysis/product/{sample}_{unit}/{sample}_{unit}.nc2.product.fna',
    output:
        merge_fna = 'analysis/product/{sample}_{unit}/{sample}_{unit}.merge.product.fna',
    shell:
        "cat {input.ext_fna}"
        " {input.nc1_fna}"
        " {input.nc2_fna}"
        " > {output.merge_fna}"

rule count_fna_product:
    input:
        merge_fna = 'analysis/product/{sample}_{unit}/{sample}_{unit}.merge.product.fna',
    output:
        clone_fna = 'analysis/logomaker/{sample}_{unit}/{sample}_{unit}.fna.clone.tsv',
    shell:
        'python {config[script_home]}/count_product.py'
        ' --ifn {input.merge_fna}'
        ' --ofn {output.clone_fna}'

rule make_fna_logomaker:
    input:
        fna = 'analysis/product/{sample}_{unit}/{sample}_{unit}.merge.product.fna' 
    output:
        fna_png = 'analysis/logomaker/{sample}_{unit}/{sample}_{unit}.fna.png',
        fna_prob = 'analysis/logomaker/{sample}_{unit}/{sample}_{unit}.fna.probability.tsv',
        fna_cnt = 'analysis/logomaker/{sample}_{unit}/{sample}_{unit}.fna.count.tsv',
        fna_fa = 'analysis/logomaker/{sample}_{unit}/{sample}_{unit}.fna'
    params:
        prefix = 'analysis/logomaker/{sample}_{unit}/{sample}_{unit}.fna'
    shell:
        'python {config[script_home]}/lm.py'
        ' --infa {input.fna}'
        ' --prefix {params.prefix}'
        ' --seqtype fna'

rule merge_faa_product:
    input:
        ext_faa = 'analysis/product/{sample}_{unit}/{sample}_{unit}.ext.product.faa',
        nc1_faa = 'analysis/product/{sample}_{unit}/{sample}_{unit}.nc1.product.faa',
        nc2_faa = 'analysis/product/{sample}_{unit}/{sample}_{unit}.nc2.product.faa',
    output:
        merge_faa = 'analysis/product/{sample}_{unit}/{sample}_{unit}.merge.product.faa' 
    shell:
        "cat {input.ext_faa}"
        " {input.nc1_faa}"
        " {input.nc2_faa}"
        " > {output.merge_faa}"

rule count_faa_product:
    input:
        merge_faa = 'analysis/product/{sample}_{unit}/{sample}_{unit}.merge.product.faa',
    output:
        clone_faa = 'analysis/logomaker/{sample}_{unit}/{sample}_{unit}.faa.clone.tsv'
    shell:
        'python {config[script_home]}/count_product.py'
        ' --ifn {input.merge_faa}'
        ' --ofn {output.clone_faa}'

rule make_faa_logomaker:
    input:
        faa = 'analysis/product/{sample}_{unit}/{sample}_{unit}.merge.product.faa' 
    output:
        faa_png = 'analysis/logomaker/{sample}_{unit}/{sample}_{unit}.faa.png',
        faa_prob = 'analysis/logomaker/{sample}_{unit}/{sample}_{unit}.faa.probability.tsv',
        faa_cnt = 'analysis/logomaker/{sample}_{unit}/{sample}_{unit}.faa.count.tsv',
        faa_fa = 'analysis/logomaker/{sample}_{unit}/{sample}_{unit}.faa'
    params:
        prefix = 'analysis/logomaker/{sample}_{unit}/{sample}_{unit}.faa'
    shell:
        'python {config[script_home]}/lm.py'
        ' --infa {input.faa}'
        ' --prefix {params.prefix}'
        ' --seqtype faa'


