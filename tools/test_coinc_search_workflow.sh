#!/bin/bash
set -e
echo -e "\\n>> [`date`] Testing pycbc_make_coinc_search_workflow"

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

echo -e "\\n>> [`date`] Building test workflow $WORKFLOWNAME"

pycbc_make_coinc_search_workflow \
--workflow-name ${WORKFLOW_NAME} --output-dir output \
--config-files \
${CONFIG_PATH}/analysis.ini \
${CONFIG_PATH}/data_O1.ini \
${CONFIG_PATH}/plotting.ini \
${CONFIG_PATH}/injections_minimal.ini \
${CONFIG_PATH}/executables.ini \
${CONFIG_PATH}/gps_times_O1_analysis_1.ini \
--config-overrides \
"results_page:output-path:../../../html"

pushd output

for workflow in *.dax
do
  echo -e "\\n>> [`date`] Validating workflow $workflow"
  pegasus-dax-validator $workflow
done

echo -e "\\n>> [`date`] Planning test workflow"
pycbc_submit_dax \
  --force-no-accounting-group \
  --dax ${WORKFLOW_NAME}.dax \
  --no-create-proxy \
  --no-submit \
  --no-grid

popd
popd

echo -e "\\n>> [`date`] Test workflow validation complete"

exit 0
