# Copyright (C) 2015  Alex Nitz
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
"""
This module is responsible for setting up plotting jobs.
https://ldas-jobs.ligo.caltech.edu/~cbc/docs/pycbc/NOTYETCREATED.html
"""
from urllib.request import pathname2url
from urllib.parse import urljoin
from pycbc.workflow.core import File, FileList, makedir, Executable

def excludestr(tags, substr):
    if substr is None:
        return tags
    if isinstance(substr, list):
        if len(substr) > 1:
            tags = excludestr(tags, substr[1:])
        substr = substr[0]
    return [tag for tag in tags if substr not in tag]

def requirestr(tags, substr):
    if substr is None:
        return tags
    return [tag for tag in tags if substr in tag]

class PlotExecutable(Executable):
    """ plot executable
    """
    current_retention_level = Executable.FINAL_RESULT

    # plots and final results should get the highest priority
    # on the job queue
    def create_node(self, **kwargs):
        node = Executable.create_node(self, **kwargs)
        node.set_priority(1000)
        return node

def make_template_plot(workflow, bank_file, out_dir,bins=None, 
                       tags=None):
    tags = [] if tags is None else tags
    makedir(out_dir)
    node = PlotExecutable(workflow.cp, 'plot_bank', ifos=workflow.ifos,
                          out_dir=out_dir, tags=tags).create_node()
    node.add_input_opt('--bank-file', bank_file)

    if workflow.cp.has_option_tags('workflow-coincidence', 'background-bins', tags=tags):
        bins = workflow.cp.get_opt_tags('workflow-coincidence', 'background-bins', tags=tags)
    if bins:
        node.add_opt('--background-bins', bins)

    node.new_output_file_opt(workflow.analysis_time, '.png', '--output-file')
    workflow += node
    return node.output_files[0]

def make_range_plot(workflow, psd_files, out_dir, exclude=None, require=None,
                   tags=None):
    tags = [] if tags is None else tags
    makedir(out_dir)
    secs = requirestr(workflow.cp.get_subsections('plot_range'), require)
    secs = excludestr(secs, exclude)
    files = FileList([])
    for tag in secs:
        node = PlotExecutable(workflow.cp, 'plot_range', ifos=workflow.ifos,
                              out_dir=out_dir, tags=[tag] + tags).create_node()
        node.add_input_list_opt('--psd-files', psd_files)
        node.new_output_file_opt(workflow.analysis_time, '.png', '--output-file')
        workflow += node
        files += node.output_files
    return files

def make_spectrum_plot(workflow, psd_files, out_dir, tags=None,
                       hdf_group=None, precalc_psd_files=None):
    tags = [] if tags is None else tags
    makedir(out_dir)
    node = PlotExecutable(workflow.cp, 'plot_spectrum', ifos=workflow.ifos,
                          out_dir=out_dir, tags=tags).create_node()
    node.add_input_list_opt('--psd-files', psd_files)
    node.new_output_file_opt(workflow.analysis_time, '.png', '--output-file')

    if hdf_group is not None:
        node.add_opt('--hdf-group', hdf_group)
    if precalc_psd_files is not None and len(precalc_psd_files) == 1:
        node.add_input_list_opt('--psd-file', precalc_psd_files)

    workflow += node
    return node.output_files[0]

def make_segments_plot(workflow, seg_files, out_dir, tags=None):
    tags = [] if tags is None else tags
    makedir(out_dir)
    node = PlotExecutable(workflow.cp, 'plot_segments', ifos=workflow.ifos,
                         out_dir=out_dir, tags=tags).create_node()
    node.add_input_list_opt('--segment-files', seg_files)
    node.new_output_file_opt(workflow.analysis_time, '.html', '--output-file')
    workflow += node

def make_gating_plot(workflow, insp_files, out_dir, tags=None):
    tags = [] if tags is None else tags
    makedir(out_dir)
    node = PlotExecutable(workflow.cp, 'plot_gating', ifos=workflow.ifos,
                          out_dir=out_dir, tags=tags).create_node()
    node.add_input_list_opt('--input-file', insp_files)
    node.new_output_file_opt(workflow.analysis_time, '.html', '--output-file')
    workflow += node

def make_throughput_plot(workflow, insp_files, out_dir, tags=None):
    tags = [] if tags is None else tags
    makedir(out_dir)
    node = PlotExecutable(workflow.cp, 'plot_throughput', ifos=workflow.ifos,
                          out_dir=out_dir, tags=tags).create_node()
    node.add_input_list_opt('--input-file', insp_files)
    node.new_output_file_opt(workflow.analysis_time, '.png', '--output-file')
    workflow += node

def make_foreground_table(workflow, trig_file, bank_file, out_dir,
                          singles=None, extension='.html', tags=None,
                          hierarchical_level=None):

    if hierarchical_level is not None and tags:
        tags = [("HIERARCHICAL_LEVEL_{:02d}".format(
                hierarchical_level))] + tags
    elif hierarchical_level is not None and not tags:
        tags = ["HIERARCHICAL_LEVEL_{:02d}".format(hierarchical_level)]
    elif hierarchical_level is None and not tags:
        tags = []

    makedir(out_dir)
    exe = PlotExecutable(workflow.cp, 'page_foreground',
                         ifos=trig_file.ifo_list,
                         out_dir=out_dir, tags=tags)
    node = exe.create_node()
    node.add_input_opt('--bank-file', bank_file)
    node.add_input_opt('--trigger-file', trig_file)
    if hierarchical_level is not None:
        node.add_opt('--use-hierarchical-level', hierarchical_level)
    if singles is not None:
        node.add_input_list_opt('--single-detector-triggers', singles)
    node.new_output_file_opt(bank_file.segment, extension, '--output-file')
    workflow += node
    return node.output_files[0]

def make_sensitivity_plot(workflow, inj_file, out_dir, exclude=None,
                         require=None, tags=None):
    tags = [] if tags is None else tags
    makedir(out_dir)
    secs = requirestr(workflow.cp.get_subsections('plot_sensitivity'), require)
    secs = excludestr(secs, exclude)
    files = FileList([])
    for tag in secs:
        node = PlotExecutable(workflow.cp, 'plot_sensitivity', ifos=workflow.ifos,
                    out_dir=out_dir, tags=[tag] + tags).create_node()
        node.add_input_opt('--injection-file', inj_file)
        node.new_output_file_opt(inj_file.segment, '.png', '--output-file')
        workflow += node
        files += node.output_files
    return files

def make_coinc_snrchi_plot(workflow, inj_file, inj_trig, stat_file, trig_file,
                          out_dir, exclude=None, require=None, tags=None):
    tags = [] if tags is None else tags
    makedir(out_dir)
    secs = requirestr(workflow.cp.get_subsections('plot_coinc_snrchi'), require)
    secs = excludestr(secs, exclude)
    files = FileList([])
    for tag in secs:
        exe = PlotExecutable(workflow.cp, 'plot_coinc_snrchi',
                             ifos=inj_trig.ifo_list,
                             out_dir=out_dir, tags=[tag] + tags)
        node = exe.create_node()
        node.add_input_opt('--found-injection-file', inj_file)
        node.add_input_opt('--single-injection-file', inj_trig)
        node.add_input_opt('--coinc-statistic-file', stat_file)
        node.add_input_opt('--single-trigger-file', trig_file)
        node.new_output_file_opt(inj_file.segment, '.png', '--output-file')
        workflow += node
        files += node.output_files
    return files

def make_inj_table(workflow, inj_file, out_dir, missed=False, singles=None,
                  tags=None):
    tags = [] if tags is None else tags
    makedir(out_dir)
    node = PlotExecutable(workflow.cp, 'page_injections', ifos=workflow.ifos,
                    out_dir=out_dir, tags=tags).create_node()

    node.add_input_opt('--injection-file', inj_file)
    if missed:
        node.add_opt('--show-missed')

    if singles is not None:
        node.add_multiifo_input_list_opt('--single-trigger-files', singles)

    node.new_output_file_opt(inj_file.segment, '.html', '--output-file')
    workflow += node
    return node.output_files[0]

def make_seg_table(workflow, seg_files, seg_names, out_dir, tags=None,
                  title_text=None, description=None):
    """ Creates a node in the workflow for writing the segment summary
    table. Returns a File instances for the output file.
    """
    seg_files = list(seg_files)
    seg_names = list(seg_names)
    if tags is None: tags = []
    makedir(out_dir)
    node = PlotExecutable(workflow.cp, 'page_segtable', ifos=workflow.ifos,
                    out_dir=out_dir, tags=tags).create_node()
    node.add_input_list_opt('--segment-files', seg_files)
    quoted_seg_names = []
    for s in seg_names:
        quoted_seg_names.append("'" + s + "'")
    node.add_opt('--segment-names', ' '.join(quoted_seg_names))
    node.add_opt('--ifos', ' '.join(workflow.ifos))
    if description:
        node.add_opt('--description', "'" + description + "'")
    if title_text:
        node.add_opt('--title-text', "'" + title_text + "'")
    node.new_output_file_opt(workflow.analysis_time, '.html', '--output-file')
    workflow += node
    return node.output_files[0]

def make_veto_table(workflow, out_dir, vetodef_file=None, tags=None):
    """ Creates a node in the workflow for writing the veto_definer
    table. Returns a File instances for the output file.
    """
    if vetodef_file is None:
        if not workflow.cp.has_option_tags("workflow-segments",
                                           "segments-veto-definer-file", []):
            return None
        vetodef_file = workflow.cp.get_opt_tags("workflow-segments",
                                           "segments-veto-definer-file", [])
        file_url = urljoin('file:', pathname2url(vetodef_file))
        vdf_file = File(workflow.ifos, 'VETO_DEFINER',
                        workflow.analysis_time, file_url=file_url)
        vdf_file.add_pfn(file_url, site='local')
    else:
        vdf_file = vetodef_file

    if tags is None: tags = []
    makedir(out_dir)
    node = PlotExecutable(workflow.cp, 'page_vetotable', ifos=workflow.ifos,
                    out_dir=out_dir, tags=tags).create_node()
    node.add_input_opt('--veto-definer-file', vdf_file)
    node.new_output_file_opt(workflow.analysis_time, '.html', '--output-file')
    workflow += node
    return node.output_files[0]

def make_seg_plot(workflow, seg_files, out_dir, seg_names=None, tags=None):
    """ Creates a node in the workflow for plotting science, and veto segments.
    """

    seg_files = list(seg_files)
    if tags is None: tags = []
    makedir(out_dir)
    node = PlotExecutable(workflow.cp, 'page_segplot', ifos=workflow.ifos,
                    out_dir=out_dir, tags=tags).create_node()
    node.add_input_list_opt('--segment-files', seg_files)
    quoted_seg_names = []
    for s in seg_names:
        quoted_seg_names.append("'" + s + "'")
    node.add_opt('--segment-names', ' '.join(quoted_seg_names))
    node.add_opt('--ifos', ' '.join(workflow.ifos))
    node.new_output_file_opt(workflow.analysis_time, '.html', '--output-file')
    workflow += node
    return node.output_files[0]

def make_ifar_plot(workflow, trigger_file, out_dir, tags=None,
                   hierarchical_level=None, executable='page_ifar'):
    """ Creates a node in the workflow for plotting cumulative histogram
    of IFAR values.
    """

    if hierarchical_level is not None and tags:
        tags = [("HIERARCHICAL_LEVEL_{:02d}".format(
                hierarchical_level))] + tags
    elif hierarchical_level is not None and not tags:
        tags = ["HIERARCHICAL_LEVEL_{:02d}".format(hierarchical_level)]
    elif hierarchical_level is None and not tags:
        tags = []

    makedir(out_dir)
    exe = PlotExecutable(workflow.cp, executable, ifos=trigger_file.ifo_list,
                         out_dir=out_dir, tags=tags)
    node = exe.create_node()
    node.add_input_opt('--trigger-file', trigger_file)
    if hierarchical_level is not None:
        node.add_opt('--use-hierarchical-level', hierarchical_level)
    node.new_output_file_opt(workflow.analysis_time, '.png', '--output-file')
    workflow += node
    return node.output_files[0]

def make_snrchi_plot(workflow, trig_files, veto_file, veto_name,
                     out_dir, exclude=None, require=None, tags=None):
    tags = [] if tags is None else tags
    makedir(out_dir)
    secs = requirestr(workflow.cp.get_subsections('plot_snrchi'), require)
    secs = excludestr(secs, exclude)
    files = FileList([])
    for tag in secs:
        for trig_file in trig_files:
            exe = PlotExecutable(workflow.cp, 'plot_snrchi',
                                 ifos=trig_file.ifo_list,
                                 out_dir=out_dir,
                                 tags=[tag] + tags)
            node = exe.create_node()

            node.set_memory(15000)
            node.add_input_opt('--trigger-file', trig_file)
            if veto_file is not None:
                node.add_input_opt('--veto-file', veto_file)
                node.add_opt('--segment-name', veto_name)
            node.new_output_file_opt(trig_file.segment, '.png', '--output-file')
            workflow += node
            files += node.output_files
    return files

def make_foundmissed_plot(workflow, inj_file, out_dir, exclude=None,
                         require=None, tags=None):
    if tags is None:
        tags = []
    makedir(out_dir)
    secs = requirestr(workflow.cp.get_subsections('plot_foundmissed'), require)
    secs = excludestr(secs, exclude)
    files = FileList([])
    for tag in secs:
        exe = PlotExecutable(workflow.cp, 'plot_foundmissed', ifos=workflow.ifos,
                    out_dir=out_dir, tags=[tag] + tags)
        node = exe.create_node()
        ext = '.html' if exe.has_opt('dynamic') else '.png'
        node.add_input_opt('--injection-file', inj_file)
        node.new_output_file_opt(inj_file.segment, ext, '--output-file')
        workflow += node
        files += node.output_files
    return files

def make_snrratehist_plot(workflow, bg_file, out_dir, closed_box=False,
                          tags=None, hierarchical_level=None):

    if hierarchical_level is not None and tags:
        tags = [("HIERARCHICAL_LEVEL_{:02d}".format(
                hierarchical_level))] + tags
    elif hierarchical_level is not None and not tags:
        tags = ["HIERARCHICAL_LEVEL_{:02d}".format(hierarchical_level)]
    elif hierarchical_level is None and not tags:
        tags = []

    makedir(out_dir)
    exe = PlotExecutable(workflow.cp, 'plot_snrratehist',
                         ifos=bg_file.ifo_list,
                         out_dir=out_dir, tags=tags)
    node = exe.create_node()
    node.add_input_opt('--trigger-file', bg_file)
    if hierarchical_level is not None:
        node.add_opt('--use-hierarchical-level', hierarchical_level)

    if closed_box:
        node.add_opt('--closed-box')

    node.new_output_file_opt(bg_file.segment, '.png', '--output-file')
    workflow += node
    return node.output_files[0]

def make_snrifar_plot(workflow, bg_file, out_dir, closed_box=False,
                     cumulative=True, tags=None, hierarchical_level=None):

    if hierarchical_level is not None and tags:
        tags = [("HIERARCHICAL_LEVEL_{:02d}".format(
                hierarchical_level))] + tags
    elif hierarchical_level is not None and not tags:
        tags = ["HIERARCHICAL_LEVEL_{:02d}".format(hierarchical_level)]
    elif hierarchical_level is None and not tags:
        tags = []

    makedir(out_dir)
    exe = PlotExecutable(workflow.cp, 'plot_snrifar', ifos=bg_file.ifo_list,
                         out_dir=out_dir, tags=tags)
    node = exe.create_node()
    node.add_input_opt('--trigger-file', bg_file)
    if hierarchical_level is not None:
        node.add_opt('--use-hierarchical-level', hierarchical_level)

    if closed_box:
        node.add_opt('--closed-box')

    if not cumulative:
        node.add_opt('--not-cumulative')

    node.new_output_file_opt(bg_file.segment, '.png', '--output-file')
    workflow += node
    return node.output_files[0]

def make_results_web_page(workflow, results_dir, template='orange',
                          explicit_dependencies=None):
    template_path = 'templates/'+template+'.html'

    out_dir = workflow.cp.get('results_page', 'output-path')
    makedir(out_dir)
    node = PlotExecutable(workflow.cp, 'results_page', ifos=workflow.ifos,
                out_dir=out_dir).create_node()
    node.add_opt('--plots-dir', results_dir)
    node.add_opt('--template-file', template_path)
    workflow += node
    if explicit_dependencies is not None:
        for dep in explicit_dependencies:
            workflow.add_explicit_dependancy(dep, node)

def make_single_hist(workflow, trig_file, veto_file, veto_name,
                     out_dir, bank_file=None, exclude=None,
                     require=None, tags=None):
    tags = [] if tags is None else tags
    makedir(out_dir)
    secs = requirestr(workflow.cp.get_subsections('plot_hist'), require)
    secs = excludestr(secs, exclude)
    files = FileList([])
    for tag in secs:
        node = PlotExecutable(workflow.cp, 'plot_hist',
                    ifos=trig_file.ifo,
                    out_dir=out_dir,
                    tags=[tag] + tags).create_node()
        if veto_file is not None:
            node.add_opt('--segment-name', veto_name)
            node.add_input_opt('--veto-file', veto_file)
        node.add_input_opt('--trigger-file', trig_file)
        if bank_file:
            node.add_input_opt('--bank-file', bank_file)
        node.new_output_file_opt(trig_file.segment, '.png', '--output-file')
        workflow += node
        files += node.output_files
    return files

def make_binned_hist(workflow, trig_file, veto_file, veto_name,
                     out_dir, bank_file, exclude=None,
                     require=None, tags=None):
    tags = [] if tags is None else tags
    makedir(out_dir)
    secs = requirestr(workflow.cp.get_subsections('plot_binnedhist'), require)
    secs = excludestr(secs, exclude)
    files = FileList([])
    for tag in secs:
        node = PlotExecutable(workflow.cp, 'plot_binnedhist',
                    ifos=trig_file.ifo,
                    out_dir=out_dir,
                    tags=[tag] + tags).create_node()
        node.add_opt('--ifo', trig_file.ifo)
        if veto_file is not None:
            node.add_opt('--veto-segment-name', veto_name)
            node.add_input_opt('--veto-file', veto_file)
        node.add_input_opt('--trigger-file', trig_file)
        node.add_input_opt('--bank-file', bank_file)
        node.new_output_file_opt(trig_file.segment, '.png', '--output-file')
        workflow += node
        files += node.output_files
    return files

def make_singles_plot(workflow, trig_files, bank_file, veto_file, veto_name,
                     out_dir, exclude=None, require=None, tags=None):
    tags = [] if tags is None else tags
    makedir(out_dir)
    secs = requirestr(workflow.cp.get_subsections('plot_singles'), require)
    secs = excludestr(secs, exclude)
    files = FileList([])
    for tag in secs:
        for trig_file in trig_files:
            node = PlotExecutable(workflow.cp, 'plot_singles',
                        ifos=trig_file.ifo,
                        out_dir=out_dir,
                        tags=[tag] + tags).create_node()

            node.set_memory(15000)
            node.add_input_opt('--bank-file', bank_file)
            if veto_file is not None:
                node.add_input_opt('--veto-file', veto_file)
                node.add_opt('--segment-name', veto_name)
            node.add_opt('--detector', trig_file.ifo)
            node.add_input_opt('--single-trig-file', trig_file)
            node.new_output_file_opt(trig_file.segment, '.png', '--output-file')
            workflow += node
            files += node.output_files
    return files

def make_dq_trigger_rate_plot(workflow, dq_files, out_dir, tags=None):
    tags = [] if tags is None else tags
    makedir(out_dir)
    files = FileList([])
    for dq_file in dq_files:
        if workflow.cp.has_option_tags('bin_trigger_rates_dq',
                                       'background-bins', tags=tags):
            background_bins = \
                          workflow.cp.get_opt_tags('bin_trigger_rates_dq',
                                              'background-bins', tags=tags)
            bin_names = [tuple(bbin.split(':'))[0] for bbin
                                             in background_bins.split(' ')]
        else: bin_names = ['all_bin']
        for bbin in bin_names:
            plot_tags = [bbin] + tags
            node = PlotExecutable(workflow.cp, 'plot_dq_likelihood_vs_time',
                        ifos=dq_file.ifo,
                        out_dir=out_dir,
                        tags=plot_tags).create_node()
            node.add_opt('--ifo', dq_file.ifo)
            node.add_opt('--background-bin', bbin)
            node.add_input_opt('--dq-file', dq_file)
            node.new_output_file_opt(dq_file.segment, '.png', '--output-file')
            workflow += node
            files += node.output_files
    return files

def make_dq_percentile_plot(workflow, dq_files, out_dir, tags=None):
    tags = [] if tags is None else tags
    makedir(out_dir)
    files = FileList([])
    for dq_file in dq_files:
        if workflow.cp.has_option_tags('bin_trigger_rates_dq',
                                       'background-bins', tags=tags):
            background_bins = \
                          workflow.cp.get_opt_tags('bin_trigger_rates_dq',
                                              'background-bins', tags=tags)
            bin_names = [tuple(bbin.split(':'))[0] for bbin
                                             in background_bins.split(' ')]
        else: bin_names = ['all_bin']
        for bbin in bin_names:
            plot_tags = [bbin] + tags
            node = PlotExecutable(workflow.cp, 'plot_dq_percentiles',
                        ifos=dq_file.ifo,
                        out_dir=out_dir,
                        tags=plot_tags).create_node()
            node.add_opt('--ifo', dq_file.ifo)
            node.add_opt('--background-bin', bbin)
            node.add_input_opt('--dq-file', dq_file)
            node.new_output_file_opt(dq_file.segment, '.png', '--output-file')
            workflow += node
            files += node.output_files
    return files
