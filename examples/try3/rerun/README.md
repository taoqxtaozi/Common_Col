## Rerun this example

Step 1. Generate the workflow plan:

```bash
bash ../../../make_plan.sh \
  --tree-file aln.tre \
  --reference GALGA \
  --paths aln.txt \
  --guide-t-threshold 0.6 \
  --threads 15 \
  --common_workers 5
```

Step 2. Execute the planned workflow.

The command above will generate an `instruction.txt` file. Follow the guidance in `instruction.txt` and run all commands listed there. Only the lines that do **not** start with `##` need to be executed.