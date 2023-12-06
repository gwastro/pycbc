""" A minimal pycbc workflow example """

import os, argparse
import pycbc
import pycbc.workflow as wf

# Standard worfklow configparser options. These give you the ability
# to supply config files on the command line, have the ability to
# override parameters, or to choose if to plan or directly submit
# the workflow
parser = argparse.ArgumentParser(description=__doc__[1:])
parser.add_argument('--verbose', action='count')
wf.add_workflow_command_line_group(parser)
wf.add_workflow_settings_cli(parser)
args = parser.parse_args()

# Set up standard logging configuration
# This is good practice for most pycbc scripts
pycbc.init_logging(args.verbose)

# This creates a file handle object to an existing file. The 
# file will typically be copied to the local directory. URLs are
# also acceptable and if provided, the file will be downloaded.
# This is useful if you have input files which come from outside
# the workflow.
input_file = wf.resolve_url_to_file(os.path.abspath('test_input.txt'))

# Create an instance of the workflow and provide it a name. Workflows
# can be contained within other workflow and each must havea a unique
# name.
workflow = wf.Workflow(args, 'cont')

##### Executable that uses arguments e.g. cp source_file dest_file
# Create an instance of an executable object. It must have a unique name.
# This name determines what keywords are needed in teh [executables] section
# of the configuration file, and any options provided in section named 
# [exe1] as in the example here will automatically be passed to all uses
# of this executable.
exe1 = wf.Executable(workflow.cp, 'argument_exe')

# The executable has a basic node creation method. The 'node' is a specific
# job you want to run using this executable. It may have options unique to the
# specific use. In this case, we have an input argument and an output argument.
node1 = exe1.create_node()

# add the input file as the first argument of the executable
node1.add_input_arg(input_file)

# Create a file handle for the future *output* of this job. Every file in
# a workflow must have a unique name. This also demonstrates making
# explicit file names and an explicit output folder. Typically,
# you may want to let the workflow decide these locations / names
out1 = wf.File.from_path(os.path.abspath("test_output.txt"))

# add the output file as the second argument of the executable
node1.add_output_arg(out1)

# Once we are finished setting up the node, we must add it to the workflow
workflow += node1

##### Executable that uses options. 
# This works the same as above, but we'll use different methods to add
# the options. We can also demonstrate the ability to automatically generate
# the file names for executable outputs.
exe2 = wf.Executable(workflow.cp, 'option_exe')
for i in range(2):
    node2 = exe2.create_node()
    
    # add the input file using a named option, note that you can
    # use as inputs files that were previously created by a another
    # job in the workflow. This will automatically have the previous job
    # run first and then pass the file it creates along to the later job.
    node2.add_input_opt('--input-file', out1)
    
    # add an option with value
    node2.add_opt('--append-count', str(i))

    # add the output file as the second argument of the executable
    # The filename must be unique and is built from the executable name
    # and any tags (set of strings) associated with the node or file.
    # We can set some which will be incorporated into the filename to make
    # sure it is unique. 
    node2.new_output_file_opt(workflow.analysis_time, '.txt',
                              "--output-file", tags=[str(i)])
    
    # Once we are finished setting up the node, we must add it to the workflow
    workflow += node2

# Save the workflow, if direction submission to your cluster was requested
# on the command line, this will also start that process automatically
workflow.save()
