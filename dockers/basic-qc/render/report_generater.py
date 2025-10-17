import io
import os
import matplotlib.pyplot as plt
import pandas as pd
import base64
from jinja2 import Environment, FileSystemLoader
import jinja2.exceptions

class ReportGenerator():
    """Generates HTML reports from report data.

    Generation takes a report data object and an HTML template name and populates
    the template.

    Formatting are methods to generate a styled HTML representation of tables
    and plots for embedding in the report.
    """

    # --- Initialization --- #
    def __init__(self, template_name, path_output):

        # assertions
        assert path_output.endswith(".html"), "Output path must end with .html"
        os.makedirs(os.path.dirname(path_output), exist_ok=True)

        # parameters
        self.template_name = template_name
        self.path_output = path_output
        self.report_data = None

        # get template (support local debugging)
        path_templates = "templates"
        if os.path.exists("/opt/templates"):
            path_templates = "/opt/templates"
        self.env = Environment(loader=FileSystemLoader(path_templates))
        try:
            self.template = self.env.get_template(self.template_name)
        except jinja2.exceptions.TemplateNotFound:
            raise ValueError(f"Template '{self.template_name}' not found in 'templates/' directory. Available templates: {self.env.list_templates()}")

    # --- HTML generation --- #
    def format_html(self, content, content_type):
        """Format HTML content for display."""
        if content_type == "text":
            return content
        elif content_type == "table":
            return self.style_table(content)
        elif content_type == "plot":
            return self.style_plot(content)
        else:
            raise ValueError(f"Unknown content type: {content_type}")

    def generate(self, report_data):
        # populate data
        self.report_data = report_data
        data = {}
        for placeholder_name, metric in self.report_data.mapped_fields.items():
            try:
                content = metric["generator"]()
                data[placeholder_name] = self.format_html(content, metric["type"])
            except Exception as e:
                content = f"<p style='color:red;'>Error generating content for '{placeholder_name}': {str(e)}</p>"
                data[placeholder_name] = self.format_html(content, "text")

        # render HTML
        html_rendered = self.template.render(**data)
        with open(self.path_output, "w") as f:
            f.write(html_rendered)

        # done
        print(f"Generated HTML report at '{self.path_output}'")

    # --- Styling methods for tables and plots --- #
    def highlight_value(row):
        """Conditional styling for DataFrame rows. If 'Warn' is True, make 'Value' red and bold. If 'Warn' is False, make 'Source & Implication' light grey."""
        styles = []
        for col in row.index:
            if col == "Value" and row["Highlight"]:
                styles.append("color:red; font-weight:bold;")
            elif col == "Source & Implication" and not row["Highlight"]:
                styles.append("color:lightgrey;")
            else:
                styles.append("") # Important: Return an empty style for other cells
        return styles

    @classmethod
    def style_table(cls, df):
        """Style the metrics table for the HTML report."""
        # format units
        df = df.copy()
        hide_columns = ["Highlight", "Unit"]
        hide_columns = [col for col in hide_columns if col in df.columns]

        # apply styling
        styled = (
            df.style
            .set_properties(**{
                "padding": "6px",
                "border": "1px solid #ccc",
                "background-color": "white",
                "color": "#333",
            })
            .set_table_styles([
                {"selector": "th", "props": [
                    ("background-color", "#f7f7f7"),
                    ("padding", "6px"),
                    ("text-align", "left"),
                    ("color", "#333")]},
                {"selector": "td.col0", "props": [("text-align", "left")]},   # Source left aligned
                {"selector": "td.col1", "props": [("text-align", "left")]}    # Implication left aligned
            ])
            .apply(cls.highlight_value, axis=1)
            .hide(axis="index")          # hide row indices
            .format("{:.1f}", subset=df.select_dtypes(include="float").columns)
            .hide(subset=hide_columns, axis="columns")
        )
        return styled.to_html()

    @classmethod
    def style_plot(cls, fig):
        """Convert a matplotlib figure to a Plotly figure for HTML embedding."""
        buf = io.BytesIO()
        fig.savefig(buf, format="png", bbox_inches="tight")
        fig.savefig(buf, format="png", dpi=200, bbox_inches="tight")
        plt.close(fig)
        buf.seek(0)
        plot_distributions_png = base64.b64encode(buf.read()).decode("ascii")
        return f"<img alt=\"Plot\" src=\"data:image/png;base64,{plot_distributions_png}\" loading=\"lazy\">"
