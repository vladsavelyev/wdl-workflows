# Configuration file templates to run WDL workflows with Cromwell

### Installation

[Cromwell](https://cromwell.readthedocs.io/) is a tool used to execute workflows written
in [WDL](https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md)
(Workflow Description Language). These workflows can be executed in a local, HPC, or
cloud (e.g. GCP) environment.

To install Cromwell, you can use conda:

```
conda create -n cromwell -c conda-forge cromwell
conda activate cromwell
```

### Configuration

Cromwell can be configured via two files passed to the `-Dconfig.file` and `--options`
command line arguments. We provide templates for both files, which are suitable for
running workflows on Google Cloud using
[Google Cloud Life Sciences](https://cromwell.readthedocs.io/en/stable/tutorials/PipelinesApi101/)
(previously known as PAPI - Pipelines API).

You must enable the Life Sciences API for your project if using the Cromwell life sciences config. You can do this from the following page (by ensuring you've selected the correct project): https://console.cloud.google.com/apis/library/lifesciences.googleapis.com.

**Important**: Make sure to schedule your VM workers in a region that's colocated with
your data buckets, to avoid incurring high network egress costs. Don't move large files
like BAMs between continents: copying 1 TB of data from the US to Australia costs 190
USD. Adjust the `default-zones`
attribute in the template if necessary.

To run Cromwell, first edit `cromwell.template.conf` to replace the following values in
angle brackets:

* `<project>`: the Google Cloud project ID (eg: project-name-12312).
* `<bucket>`: a Google Cloud Storage bucket name to store executions.
* (optional) `<mysql-password>`: MySQL password for a locally running MySQL server that
  will be used to track Cromwell executions in order to allow restarts of incomplete
  runs. You can comment out the entire `database` section if you don't need that
  functionality.

Also edit `options.template.json`:

* `<bucket>`: a Google Cloud Storage bucket name to store outputs of successfully
  finished executions
  
#### Authentication

By default, Cromwell will use the account that's currently authenticated with `gsutil` on the current computer. It's possible to configure Cromwell to use a service account, see the [Cromwell: Google Backend](https://cromwell.readthedocs.io/en/stable/backends/Google/) for more information about configuring Cromwell to authenticate using a service account.
  
### Running workflows

Finally, to run a workflow `workflow.wdl` with inputs in `inputs.json`, use the
following command (assuming the edited configuration files are saved
under `cromwell.conf` and `options.json`):

```
cromwell -Dconfig.file=cromwell.conf run workflow.wdl --inputs inputs.json --options options.json &
```

Make sure to keep `&` in the end to run the process in the background, otherwise you
might accidentally interrupt the execution. Use `fg` to bring the process back to the
foreground.



## WARP inputs

The `warp-input-templates/` folder contains templates that can be modified to use with
the [germline variant calling WARP workflows](https://github.com/populationgenomics/warp/blob/master/pipelines/broad/dna_seq/germline/).
You'll have to replace the input parameters at the top. Specifically:

* `<sample-name>` and `<bam-location>` for the single-sample workflows `WGSFromBam`
  and `ExomeFromBam`,
* parameters in the `sample_and_fastqs` section for `WGSFromFastq`,
* or the pointer to a `sample_map` file location for `ExomeMultipleSamplesFromBam`
  , `WGSMultipleSamplesFromBam` or `WGSMultipleSamplesFromFastq`, where the sample map
  is a tab-separated file with 2 columns: the sample name, and the input file location.
  For example:

```bash
NA11843	gs://genomics-public-data/ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/NA11843/alignment/NA11843.mapped.ILLUMINA.bwa.CEU.low_coverage.20120522.bam
```

All the reference-data inputs are pre-filled to point to the Broad public genomics
buckets.


## Examples

To run
a [WGS single-sample germline variant calling WDL workflow](https://github.com/populationgenomics/warp/blob/start_from_mapped_bam/pipelines/broad/dna_seq/germline/single_sample/) (
which is based on [Broad WARP](https://github.com/broadinstitute/warp/)) on one sample
using the data from `gs://genomics-public-data`, run the following commands:

```
git clone https://github.com/populationgenomics/fewgenomes
git clone https://github.com/populationgenomics/warp
SAMPLE=HG00272
cromwell -Dconfig.file=cromwell.conf run \
    warp/pipelines/broad/dna_seq/germline/single_sample/wgs/WGSFromBam.wdl \ 
    --inputs fewgenomes/datasets/2genomes/wgs_bam/$SAMPLE.json \
    --options options.json
```

To run the WGS workflow on multiple samples in parallel for the entire dataset 
`6genomes`, use:

```
git clone https://github.com/populationgenomics/fewgenomes
git clone https://github.com/populationgenomics/warp
cromwell -Dconfig.file=cromwell.conf run \
    warp/pipelines/broad/dna_seq/germline/single_sample/wgs/WGSMultipleSamplesFromBam.wdl \ 
    --inputs fewgenomes/datasets/6genomes/6genomes-wgs_bam.json \
    --options options.json
```

Also
see [populationgenomics/fewgenomes](https://github.com/populationgenomics/fewgenomes)
for more details on how the input JSON was generated.
