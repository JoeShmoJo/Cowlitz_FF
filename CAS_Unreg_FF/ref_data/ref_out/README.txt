Small, representative sample output for this project -- e.g. the
result of running your pipeline against ref_in/, kept here so expected
output can be checked/diffed without regenerating the full dataset.

Keep this folder to a handful of files. Clean out test-run leftovers
before committing -- e.g. in PowerShell:
  Get-ChildItem -Path ref_data\ref_out -Exclude README.txt | Remove-Item -Recurse -Force
