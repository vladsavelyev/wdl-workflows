# Configuration file templates to run WDL workflows with Cromwell

[Cromwell](https://cromwell.readthedocs.io/) is a tool used to execute
workflows written in [WDL](https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md)
(Workflow Description Language). These workflows can be executed in a local,
HPC, or cloud (e.g. GCP) environment.

To install Cromwell, you can use conda:

```
create -n cromwell -c conda-forge cromwell
conda activate cromwell
```

Cromwell can be configured via two files passed to the `-Dconfig.file` and `--options`
command line arguments. We provide templates for both files, which are suitable
for running workflows on Google Cloud using
[Google Cloud Life Sciences](https://cromwell.readthedocs.io/en/stable/tutorials/PipelinesApi101/)
(previously known as PAPI - Pipelines API).

**Important**: Make sure to schedule your VM workers in a region that's
colocated with your data buckets, to avoid incurring high network egress
costs. Don't move large files like BAMs between continents: copying 1 TB of
data from the US to Australia costs 190 USD. Adjust the `default-zones`
attribute in the template if necessary.

To run Cromwell, first edit `cromwell.template.conf` to replace the following values in angle brackets:

* `<project>`: the Google Cloud project name.
* `<bucket>`: a Google Cloud Storage bucket name to store executions.
* (optional) `<mysql-password>`: MySQL password for a locally running MySQL server that will be used to track Cromwell executions in order to allow restarts of incomplete runs. You can comment out the entire `database` section if you don't need that functionality.

Also edit `options.template.json`:

* `<bucket>`: a Google Cloud Storage bucket name to store outputs of successfully finished executions

Finally, to run a workflow `workflow.wdl` with inputs in `inputs.json`, use the following command (assuming the edited configuration files are saved under `cromwell.conf` and `options.json`):

```
cromwell -Dconfig.file=cromwell.conf run workflow.wdl --inputs inputs.json --options options.json &
```

Make sure to keep `&` in the end to run the process in the background,
otherwise you might accidentally interrupt the execution.
Use `fg` to bring the process back to the foreground.

## WARP inputs

The `warp-inputs/` folder contains templates that can be modified to use with [WGS or WES germline variant calling WARP workflows](https://github.com/populationgenomics/warp/blob/master/pipelines/broad/dna_seq/germline/single_sample/). You'll have to replace the input parameters at the top (`<sample-name>` and `<bam-location>` for the workflows that start from BAM, or the parameters in the `sample_and_fastqs` section for the workflows that start from FASTQs); all the reference-data inputs are pre-filled to point to the Broad public genomics buckets.

## Examples

To run a [WES germline variant calling WDL workflow](https://github.com/populationgenomics/warp/blob/start_from_mapped_bam/pipelines/broad/dna_seq/germline/single_sample/) (which is based on [Broad WARP](https://github.com/broadinstitute/warp/)) on a sample NA12878 using the data from `gs://genomics-public-data`, run the following commands:

```
git clone https://github.com/populationgenomics/fewgenomes
git clone https://github.com/populationgenomics/warp
SAMPLE=NA12878
cromwell -Dconfig.file=cromwell.conf run \
    warp/pipelines/broad/dna_seq/germline/single_sample/exome/ExomeFromBam.wdl \ 
    --inputs fewgenomes/datasets/toy/exome_bam/$SAMPLE.json \
    --options options.json
```

Also see [populationgenomics/cpg-fewgenomes](https://github.com/populationgenomics/cpg-fewgenomes) for more details on how the input JSON was generated.
