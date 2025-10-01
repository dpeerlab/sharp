import os
import click
from report_generater import ReportGenerator
from report_data import ReportData

@click.group()
def render_sharp_report():
    pass

@render_sharp_report.command(no_args_is_help=True)
@click.argument("template-name")
@click.option("--sample-name", help="Name of the sample.", default="test_hashtag")
@click.option("--path-h5ad", help="Path to the input h5ad file.", default="data/hashtag/adata.h5ad")
@click.option("--path-report", help="Path to the run report YAML file.", default="data/hashtag/run_report.yaml")
@click.option("--path-reads", help="Path to the directory with read files.", default= "data/hashtag/reads")
@click.option("--path-output", help="Path to the output HTML file.", default="outputs/hashtag/hashtag_report.html")
def render(template_name, sample_name, path_h5ad, path_report, path_reads, path_output):
    # assertions
    assert os.path.exists(path_h5ad), f"Input h5ad file '{path_h5ad}' does not exist."
    assert os.path.exists(path_report), f"Run report file '{path_report}' does not exist."
    assert os.path.exists(path_reads), f"Reads directory '{path_reads}' does not exist."

    # setup report
    html_generator = ReportGenerator(
        template_name=template_name,
        path_output=path_output,
    )

    # get data
    report_data = ReportData(
        sample_name=sample_name,
        path_h5ad=path_h5ad,
        path_report=path_report,
        path_reads=path_reads
    )

    # generate report
    html_generator.generate(
        report_data=report_data
    )

if __name__ == "__main__":
    render_sharp_report()