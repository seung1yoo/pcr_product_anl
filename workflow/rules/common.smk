from snakemake.utils import validate
from snakemake.utils import Paramspace
import pandas as pd

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"

##### load config and sample sheets #####

#configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")
print(config)

samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"sample_name": str})
    .set_index("sample_name", drop=False)
    .sort_index()
)
print(samples)
validate(samples, schema="../schemas/samples.schema.yaml")

units = (
    pd.read_csv(config["units"], sep="\t", dtype={"sample_name": str, "unit_name": str})
    .set_index(["sample_name", "unit_name"], drop=False)
    .sort_index()
)
print(units)
validate(units, schema="../schemas/units.schema.yaml")

#paramspace = Paramspace(pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False))

def get_outputs():
    outputs = list()
    for sample, unit in units.index:
        outputs.append(f'analysis/flash/{sample}_{unit}/{sample}_{unit}.extendedFrags.fastq.gz')
        outputs.append(f'analysis/flash/{sample}_{unit}/{sample}_{unit}.notCombined_1.fastq.gz')
        outputs.append(f'analysis/flash/{sample}_{unit}/{sample}_{unit}.notCombined_2.fastq.gz')
        outputs.append(f'analysis/flash/{sample}_{unit}/{sample}_{unit}.hist')
        outputs.append(f'analysis/flash/{sample}_{unit}/{sample}_{unit}.histogram')

        outputs.append(f'analysis/product/{sample}_{unit}/{sample}_{unit}.ext.product.fna')
        outputs.append(f'analysis/product/{sample}_{unit}/{sample}_{unit}.ext.product.faa')
        outputs.append(f'analysis/product/{sample}_{unit}/{sample}_{unit}.ext.product.tsv')
        outputs.append(f'analysis/product/{sample}_{unit}/{sample}_{unit}.nc1.product.fna')
        outputs.append(f'analysis/product/{sample}_{unit}/{sample}_{unit}.nc1.product.faa')
        outputs.append(f'analysis/product/{sample}_{unit}/{sample}_{unit}.nc1.product.tsv')
        outputs.append(f'analysis/product/{sample}_{unit}/{sample}_{unit}.nc2.product.fna')
        outputs.append(f'analysis/product/{sample}_{unit}/{sample}_{unit}.nc2.product.faa')
        outputs.append(f'analysis/product/{sample}_{unit}/{sample}_{unit}.nc2.product.tsv')
        outputs.append(f'analysis/product/{sample}_{unit}/{sample}_{unit}.merge.product.fna')
        outputs.append(f'analysis/product/{sample}_{unit}/{sample}_{unit}.merge.product.faa')

        outputs.append(f'analysis/logomaker/{sample}_{unit}/{sample}_{unit}.fna')
        outputs.append(f'analysis/logomaker/{sample}_{unit}/{sample}_{unit}.fna.png')
        outputs.append(f'analysis/logomaker/{sample}_{unit}/{sample}_{unit}.fna.probability.tsv')
        outputs.append(f'analysis/logomaker/{sample}_{unit}/{sample}_{unit}.fna.count.tsv')
        outputs.append(f'analysis/logomaker/{sample}_{unit}/{sample}_{unit}.fna.clone.tsv')

        outputs.append(f'analysis/logomaker/{sample}_{unit}/{sample}_{unit}.faa')
        outputs.append(f'analysis/logomaker/{sample}_{unit}/{sample}_{unit}.faa.png')
        outputs.append(f'analysis/logomaker/{sample}_{unit}/{sample}_{unit}.faa.probability.tsv')
        outputs.append(f'analysis/logomaker/{sample}_{unit}/{sample}_{unit}.faa.count.tsv')
        outputs.append(f'analysis/logomaker/{sample}_{unit}/{sample}_{unit}.faa.clone.tsv')
    return(outputs)

    

    
    



    




