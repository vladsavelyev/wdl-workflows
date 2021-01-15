# Configuration file templates to run WDL workflows with Cromwell

To run cromwell, edit `cromwell.template.conf`. Replace the following values:

* `<project>` Google Cloud project name
* `<bucket>`: bucket name for executions
* `<mysql-password>`: MySQL password for a local MySQL server that is used to track Cromwell executions to allow restarts of incomplete runs. You can comment out the entire `database` section if you don't need that functionality.

Also edit `options.template.json`

* `<bucket>`: bucket name for executions and outputs
* `<local-dir-for-cromwell-logs>`: local directory to write Cromwell logs

If you don't have `cromwell` in your PATH, install a conda environment with `cromwell`:

```
create -n cromwell cromwell
conda activate cromwell
```

To run a workflow `workflow.wdl`:

```
cromwell -Dconfig.file=cromwell.conf run workflow.wdl --inputs inputs.json --options options.json &
```

Make sure to keep `&` in the end to run the process in the background, otherwise you might accidentally interrupt the execution. To bring the process to the foreground, so you could interrupt it, run `fg`.
