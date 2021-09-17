# WDL workflows and Cromwell configuration files

## WDL workflows

[WDL](https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md) (Workflow Description Language) is a language developed by the Broad Institute that allows to write genomics workflows specification. WDL describes commands to invoke tools, required computation resources, and how tools are piped together; but it abstracts from particular execution environment implementations.

This repository provides workflows for alignment and germline variant calling of whole genome and exome DNA sequencing data. The `workflows` folder contains two WDL files: `SingleSample.wdl` and `MultipleSamples.wdl`. The former one is based on the [WARP WholeGenomeGermlineSingleSample](https://github.com/broadinstitute/warp/tree/develop/pipelines/broad/dna_seq/germline/single_sample/wgs) workflow, modified to make use of [Bazam](https://github.com/ssadedin/bazam) and [biobambam2](https://github.com/gt1/biobambam2) to stream the enitre alignment process from input BAM/CRAM down to re-aligned CRAM in just single command. Additionally, the exome and WGS functionalities are combined into a single workflow: the former is chosen if `target_interval_list` and `bait_interval_list` are specified, whereas the latter is triggered if `wgs_coverage_interval_list` is specified.

The `test-inputs` folder contains examples of input JSONs for the `SingleSample.wdl` and `MultipleSamples.wdl` workflows. Edit the `SingleSample.inp` section to set your own input files (which can be a BAM with an index, a CRAM with an index, or a pair of FASTQs). For `MultipleSamples`, edit a `*.sample-map.tsv` file, and point to it in the `MultipleSamples.sample_map` section of the JSON.

All the reference-data inputs are pre-filled to point to the Broad public genomics buckets.

## Cromwell

[Cromwell](https://cromwell.readthedocs.io/) is a tool used to execute workflows written in WDL. These workflows can be executed in a local, HPC, or cloud (e.g. GCP) environment, however GCP is preferable for CPG.

### Installation

To install Cromwell, you can use conda:

```bash
create -n cromwell -c conda-forge cromwell
conda activate cromwell
```

### Prerequisites

You must enable the Life Sciences API for your project when using the provided template. You can do this from the following page; make sure you've selected the correct project: https://console.cloud.google.com/apis/library/lifesciences.googleapis.com.

Cromwell can use a MySQL database to keep track of workflow executions. This is required to enable call-caching, which lets you resume workflows by reusing intermediate results. While not strictly necessary, this option is highly recommended, partly because the alternative in-memory database easily runs out of space.

```bash
brew install mysql
```

Start a local MySQL server and connect to it:

```bash
mysql.server start
mysql -u root
```

Create a database called `cromwell`, by executing the following at the MySQL prompt:

```sql
CREATE DATABASE cromwell;
exit
```

### Configuration

Cromwell can be configured via two files passed to the `-Dconfig.file` and `--options` command line arguments. We provide templates for both files, which are suitable for running workflows on Google Cloud using [Google Cloud Life Sciences](https://cromwell.readthedocs.io/en/stable/tutorials/PipelinesApi101/) (previously known as PAPI - Pipelines API).

**Important**: Make sure to schedule your VM workers in a region that's colocated with your data buckets, to avoid incurring high network egress costs. Don't move large files like BAMs between continents: copying 1 TB of data from the US to Australia costs 190 USD. Adjust the `default-zones` attribute in the template if necessary.

To run Cromwell, first edit `cromwell.conf` to replace the following hardcoded values:

* `project = "fewgenomes"` (defined in `engine` and `backend/providers` sections): Google Cloud project ID (e.g. `project-name-12312`),
* `root = "gs://cpg-fewgenomes-test-tmp/cromwell/executions"` (defined in the `backend/providers` section): a Google Cloud Storage bucket location to store executions,
* `password = "12345678"` (defined in the `database` section): MySQL password for a locally running MySQL server that will be used to track Cromwell executions in order to allow restarts of incomplete runs. You can comment out the entire `database` section if you don't need that functionality.

Also edit `options.json` or `options-cache.json`, and replace the `gs://cpg-fewgenomes-main` bucket name with your bucket to store logs and outputs of successfully finished executions.

### Authentication

By default, Cromwell will use the account that's currently authenticated with `gcloud` on the current computer. It's possible to configure Cromwell to use a service account, see the [Cromwell: Google Backend](https://cromwell.readthedocs.io/en/stable/backends/Google/) for more information about configuring Cromwell to authenticate using a service account.


## Running workflows

Finally, to run a workflow (for example, the single sample alignment and germline variant callign workflow `workflows/SingleSample.wdl` with inputs `test-inputs/SingleSample.wgs-cram.json`), use the following command:

```bash
cromwell -Dconfig.file=cromwell.conf run workflows/SingleSample.wdl --inputs test-inputs/SingleSample.wgs-cram.json --options options.json &
```

Make sure to keep `&` in the end to run the process in the background, otherwise you might accidentally interrupt the execution. Use `fg` to bring the process back to the foreground.

### Examples

Run the single-sample exome alignment and variant calling workflow from an input BAM:

```bash
cromwell -Dconfig.file=cromwell.conf run workflows/SingleSample.wdl --inputs test-inputs/SingleSample.exome-bam.json --options options-cache.json &
```

Run multiple single-sample workflows in parallel from input pairs of FASTQs:

```bash
cromwell -Dconfig.file=cromwell.conf run workflows/MultipleSamples.wdl --inputs test-inputs/MultipleSamples.wgs-fastq.json --options options-cache.json &
```
