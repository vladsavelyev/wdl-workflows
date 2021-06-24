# WDL workflows and Cromwell configuration files

[Cromwell](https://cromwell.readthedocs.io/) is a tool used to execute workflows written
in [WDL](https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md)
(Workflow Description Language). These workflows can be executed in a local, HPC, or
cloud (e.g. GCP) environment.

To install Cromwell, you can use conda:

```
create -n cromwell -c conda-forge cromwell
conda activate cromwell
```

Cromwell can be configured via two files passed to the `-Dconfig.file` and `--options`
command line arguments. We provide templates for both files, which are suitable for
running workflows on Google Cloud using
[Google Cloud Life Sciences](https://cromwell.readthedocs.io/en/stable/tutorials/PipelinesApi101/)
(previously known as PAPI - Pipelines API).

**Important**: Make sure to schedule your VM workers in a region that's colocated with
your data buckets, to avoid incurring high network egress costs. Don't move large files
like BAMs between continents: copying 1 TB of data from the US to Australia costs 190
USD. Adjust the `default-zones`
attribute in the template if necessary.

To run Cromwell, first edit `cromwell.conf` to replace the following hardcoded values:

* `project = "fewgenomes"` (defined in `engine` and `backend/providers` sections): your Google Cloud project name,
* `root = "gs://cpg-fewgenomes-test-tmp/cromwell/executions"` (defined in the `backend/providers` section): a Google Cloud Storage bucket location to store executions,
* (optional) `password = "12345678"` (defined in the `database` section): MySQL password for a locally running MySQL server that will be used to track Cromwell executions in order to allow restarts of incomplete runs. You can comment out the entire `database` section if you don't need that  functionality.

Also edit `options.json` or `options-cache.json`, and replace the `gs://cpg-fewgenomes-main` bucket name with your bucket to store logs and outputs of successfully finished executions.

Finally, to run a workflow (for example, the single sample alignment and germline variant callign workflow `workflows/SingleSample.wdl` with inputs `test-inputs/SingleSample.wgs-cram.json`), use the following command:

```
cromwell -Dconfig.file=cromwell.conf run workflows/SingleSample.wdl --inputs test-inputs/SingleSample.wgs-cram.json --options options.json &
```

Make sure to keep `&` in the end to run the process in the background, otherwise you
might accidentally interrupt the execution. Use `fg` to bring the process back to the
foreground.

## WDL workflows

The `workflows` folder contains two WDL workflows for alignment and germline variant calling: `SingleSample.wdl` and `MultipleSamples.wdl`. The former is based on [WARP workflows]((https://github.com/populationgenomics/warp/blob/master/pipelines/broad/dna_seq/germline/) by the Broad Institute, however it makes use of Bazam and Biobamba to stream the alignment process from input BAM/CRAM into output CRAM in one command. Also, the exome and WGS functionalities are combined into a single workflow, and triggered depending whether `wgs_coverage_interval_list` is specified, or `target_interval_list` and `bait_interval_list`.

The `test-inputs` contains examples of input JSONs for the `SingleSample.wdl` and `MultipleSamples.wdl` workflows. Edit the `SingleSample.inp` section to set your own inputs. For `MultipleSamples`, also edit a `*.sample-map.tsv` file, and point to it in the `MultipleSamples.sample_map` section of the JSON.

All the reference-data inputs are pre-filled to point to the Broad public genomics
buckets.

## Examples

Run the single-sample exome alignment and variant calling workflow from an input BAM:

```bash
cromwell -Dconfig.file=cromwell.conf run workflows/SingleSample.wdl --inputs test-inputs/SingleSample.exome-bam.json --options options-cache.json
```

Run multiple single-sample workflows in parallel from input pairs of FASTQs:

```bash
cromwell -Dconfig.file=cromwell.conf run workflows/MultipleSamples.wdl --inputs test-inputs/MultipleSamples.wgs-fastq.json --options options-cache.json```
```
