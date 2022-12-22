# Primary CLIP Analysis

## Quick-start

To test the pipeline, use the associated config file and run it with the
profiles `test` and the container engine you wish to use eg. `docker`. For example:

```bash
nextflow run main.nf -profile test,docker
```

The `test` profile will auto-add all params from the `assets` folder.
