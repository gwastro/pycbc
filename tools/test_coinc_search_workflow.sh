#!/bin/bash
set -e
echo -e "\\n>> [`date`] Testing pycbc_make_offline_search_workflow"

VENV_PATH=${1}
TRAVIS_TAG=${2}

if [ "x${VENV_PATH}" == "x" ] ; then
  echo -e "\\n>> [`date`] Error: VENV_PATH was not passed to script or is empty"
  exit 1
fi

if [ "x${TRAVIS_TAG}" == "x" ] ; then
  echo -e "\\n>> [`date`] Error: TRAVIS_TAG was not passed to script or is empty"
  exit 1
fi

echo -e "\\n>> [`date`] Entering virtual environment $VENV_PATH"
source ${VENV_PATH}/bin/activate

CONFIG_PATH="https://raw.githubusercontent.com/ligo-cbc/pycbc-config/${TRAVIS_TAG}/test"
echo -e "\\n>> [`date`] Using config files from ${CONFIG_PATH}"

echo -e "\\n>> [`date`] Creating test workflow"
UUID=`uuidgen`
WORKFLOW_NAME=test-workflow-$UUID
OUTPUT_PATH=`pwd`/public_html/test_workflow/${WORKFLOW_NAME}
export LIGO_DATAFIND_SERVER='128.230.190.43:80'

mkdir $WORKFLOW_NAME
pushd $WORKFLOW_NAME

# Doesn't need to be a valid bank file, just needs to exist
echo "DUMMY BANK FILE" > bank.hdf
echo "DUMMY STAT FILE" > statHL.hdf
echo "DUMMY STAT FILE" > statLV.hdf
echo "DUMMY STAT FILE" > statHV.hdf
echo "DUMMY STAT FILE" > statHLV.hdf

echo -e "\\n>> [`date`] Building test workflow $WORKFLOWNAME"

pycbc_make_offline_search_workflow \
--workflow-name ${WORKFLOW_NAME} --output-dir output \
--config-files \
/pycbc/examples/search/analysis.ini \
/pycbc/examples/search/plotting.ini \
/pycbc/examples/search/injections_minimal.ini \
/pycbc/examples/search/executables.ini \
--config-overrides \
"results_page:output-path:../../../html"

pushd output

for workflow in *.dax
do
  echo -e "\\n>> [`date`] Validating workflow $workflow"
  pegasus-dax-validator $workflow
done

echo -e "\\n>> [`date`] Submitting test workflow"
bash start

popd
popd

echo -e "\\n>> [`date`] Test workflow validation complete"

exit 0
