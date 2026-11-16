# ruff: noqa: I001
# NOTE: pygrb_plotting_utils / pygrb_postprocessing_utils MUST be imported
# last. They pull in pycbc.io.gracedb, which does
# `from pycbc.results import generate_asd_plot, generate_snr_plot,
# source_color` -- a circular import back into this partially-initialized
# module. That only works if color/psd/snr (which define those names) have
# already been imported above. A "safe" isort/ruff I001 autofix alphabetized
# these imports once already, moving the pygrb_* imports before snr and
# breaking this exact thing -- do not let that happen again.
from pycbc.results.color import (
    ifo_color,
    source_color,
)
from pycbc.results.dq import (
    redirect_javascript,
    search_form_string,
    data_h1_string,
    data_l1_string,
    get_summary_page_link,
)
from pycbc.results.layout import (
    two_column_layout,
    single_layout,
    grouper,
    group_layout,
    SectionNumber,
)
from pycbc.results.metadata import (
    escape_table,
    unescape_table,
    html_escape,
    MetaParser,
    save_html_with_metadata,
    load_html_metadata,
    save_png_with_metadata,
    save_pdf_with_metadata,
    load_png_metadata,
    save_fig_with_metadata,
    load_metadata_from_file,
)
from pycbc.results.plot import (
    hist_overflow,
    add_style_opt_to_parser,
    set_style_from_cli,
)
from pycbc.results.psd import (
    generate_asd_plot,
)
from pycbc.results.snr import (
    generate_snr_plot,
)
from pycbc.results.str_utils import (
    drop_trailing_zeros,
    format_value,
    get_signum,
    mathjax_html_header,
)
from pycbc.results.table_utils import (
    google_table_template,
    html_table,
    static_table_template,
    static_table,
)
from pycbc.results.versioning import (
    get_library_version_info,
    get_code_version_numbers,
)
from pycbc.results.pygrb_plotting_utils import (
    contour_plotter,
    make_grb_segments_plot,
    axis_max_value,
    axis_min_value,
    pygrb_plotter,
)
from pycbc.results.pygrb_postprocessing_utils import (
    pygrb_initialize_plot_parser,
    pygrb_add_slide_opts,
    slide_opts_helper,
    pygrb_add_injmc_opts,
    pygrb_add_bestnr_opts,
    pygrb_add_null_snr_opts,
    pygrb_add_bestnr_cut_opt,
    load_data,
    apply_vetoes_to_found_injs,
    get_antenna_responses,
    get_antenna_dist_factor,
    sort_trigs,
    extract_trig_properties,
    extract_ifos,
    load_time_slides,
    load_segment_dict,
    construct_trials,
    sort_stat,
    max_median_stat,
    mc_cal_wf_errs,
    get_coinc_snr,
    template_hash_to_id,
)

__all__ = [
    "MetaParser",
    "SectionNumber",
    "add_style_opt_to_parser",
    "apply_vetoes_to_found_injs",
    "axis_max_value",
    "axis_min_value",
    "construct_trials",
    "contour_plotter",
    "data_h1_string",
    "data_l1_string",
    "drop_trailing_zeros",
    "escape_table",
    "extract_ifos",
    "extract_trig_properties",
    "format_value",
    "generate_asd_plot",
    "generate_snr_plot",
    "get_antenna_dist_factor",
    "get_antenna_responses",
    "get_code_version_numbers",
    "get_coinc_snr",
    "get_library_version_info",
    "get_signum",
    "get_summary_page_link",
    "google_table_template",
    "group_layout",
    "grouper",
    "hist_overflow",
    "html_escape",
    "html_table",
    "ifo_color",
    "load_data",
    "load_html_metadata",
    "load_metadata_from_file",
    "load_png_metadata",
    "load_segment_dict",
    "load_time_slides",
    "make_grb_segments_plot",
    "mathjax_html_header",
    "max_median_stat",
    "mc_cal_wf_errs",
    "pygrb_add_bestnr_cut_opt",
    "pygrb_add_bestnr_opts",
    "pygrb_add_injmc_opts",
    "pygrb_add_null_snr_opts",
    "pygrb_add_slide_opts",
    "pygrb_initialize_plot_parser",
    "pygrb_plotter",
    "redirect_javascript",
    "save_fig_with_metadata",
    "save_html_with_metadata",
    "save_pdf_with_metadata",
    "save_png_with_metadata",
    "search_form_string",
    "set_style_from_cli",
    "single_layout",
    "slide_opts_helper",
    "sort_stat",
    "sort_trigs",
    "source_color",
    "static_table",
    "static_table_template",
    "template_hash_to_id",
    "two_column_layout",
    "unescape_table",
]
