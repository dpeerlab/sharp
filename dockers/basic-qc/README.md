# sharp-basic-qc


## Code Overview

Reports consist of 3 main components:

1. HTML templates
2. Data extraction script
3. CLI to generate the report

Important files:

- `templates/citeseq_report.html` and `templates/hashtag_report.html` define the HTML templates for generating the reports.
- `render/report_data.py` contains the script to extract the data:
    - `__init__(...)`: The field `self.mapped_fields` renders plots, tables etc. in the respective `html` files.
    - `generate_table_metrics(...)`: Generates summary metrics tables.
    - `generate_table_warnings(...)`: Same as before, but only show if `highlight` is `True`.
- `cli/cli.py` is the command-line interface to generate the report.

## Test

### Local

Ensure that `hto` is installed

```
source run-local.sh
```

### Using `hto` docker

```
