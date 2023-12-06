# In this example, we run the workflow generation script and immediately
# submit it to the compute cluster

# The env is modified only only so that the executable in this folder 
# 'option_exe' can be used to demonstrate the config file 'which' usage.
export PATH=$PATH:$PWD

python simple.py \
--workflow-name test \
--config-files simple.ini \
--submit-now
